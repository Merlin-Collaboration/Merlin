/////////////////////////////////////////////////////////////////////////
// class KSEnergyAdjustment implementation
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/19 10:19:06 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "KSEnergyAdjustment.h"
#include "ILCDFS_IO.h"

using namespace std;

KSEnergyAdjustment::KSEnergyAdjustment(double dEr) 
: EnergyAdjustmentPolicy(), delta(dEr)
{}

void KSEnergyAdjustment::Initialise()
{
	dfs_trace(dfs_trace::level_1)<<"Initialising Klystron Shunting with "<<100*delta<<"%"<<endl;
	// get initial beam energy
	energy0 = beamrefs[0]->GetReferenceMomentum();
}

void KSEnergyAdjustment::SetActiveBeamlineSegment(DFS_Segment &seg)
{
	// Identify last klystron upstream of segment
	double energy = energy0;
	size_t nk1,nk2;
	for(size_t i=0; i<theKlystrons.size(); i++) {
		std::vector<size_t> ki;
		theKlystrons[i]->GetBeamlineIndecies(ki);
		if(ki.back() >= seg.first) {
			nk2 = i==0 ? 0 : i-1;
			break;
		}
		energy += (*theKlystrons[i]).GetVoltagePhasor().real();
	}

	// Indentify klystron range to modify
	double v = delta*energy;
	nk1=nk2+1;
	if(nk2!=0) {
		while(v>0 && nk1>0) {
			nk1--;
			v-=(*theKlystrons[nk1]).GetVoltagePhasor().real();
		}
		dEbeam = v>0 ? -v : 0; // adjust initial beam if necessary
		cKlysRange.first = nk1;
		cKlysRange.second = nk2;
	}
	else { // special case: no klystrons can be used
		dEbeam = -v;
		cKlysRange.first = 0;
		cKlysRange.second = 0;
	}
	if(dfs_trace::verbosity>=dfs_trace::level_2) {
		if(nk2==0)
			dfs_trace(dfs_trace::level_2)<<"zero klystrons modified; ";
		else {
			dfs_trace(dfs_trace::level_2)<<nk2-nk1+1<<" klystrons ("<<nk1<<'-'<<nk2<<") ";
			dfs_trace(dfs_trace::level_2)<<delta*energy-v<<" GeV; ";
		}
		dfs_trace(dfs_trace::level_2)<<"beam = "<<dEbeam<<" GeV"<<endl;
	}
};

void KSEnergyAdjustment::SetEnergyState(size_t nes)
{
	if(nes==0) {
		beamrefs[1]->SetReferenceMomentum(energy0);
		RestoreKlystrons();
	}
	else {
		if(cKlysRange.first!=0) {
			for(size_t k=cKlysRange.first; k<=cKlysRange.second; k++) {
				theKlystrons[k]->SetVoltage(0.0);
			}
		}
		beamrefs[1]->SetReferenceMomentum(energy0+dEbeam);
	}
}



