/////////////////////////////////////////////////////////////////////////
// class ConstantGradientAdjustment implementation
//
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "ConstantGradientAdjustment.h"
#include "ILCDFS_IO.h"

using namespace std;

void ConstantGradientAdjustment::Initialise()
{
	// Here we adjust the initial energy of the beam states

	dfs_trace(dfs_trace::level_1)<<"\n\nInitialising initial beam energy states\n"<<endl;

	for(size_t n=0; n<estates.size(); n++) {
		double p0 = beamrefs[n]->GetReferenceMomentum();
		double p=p0*(1+estates[n].dEbeam);
		beamrefs[n]->SetReferenceMomentum(p);
		dfs_trace(dfs_trace::level_2)<<"state "<<n<<": p = "<<p<<" GeV/c"<<endl;
	}
}

void ConstantGradientAdjustment::AddEnergyState(double dEk, double dEb)
{
	dfs_trace(dfs_trace::level_3)<<"Adding energy state ("<<dEk<<", "<<dEb<<")"<<endl;
	estates.push_back(EnergyState(dEk,dEb));
}

void ConstantGradientAdjustment::SetEnergyState(size_t nes)
{
	const EnergyState& es = estates[nes];
	// Note that the non-zero initial beam energy states
	// have already been set during the call to Initialise().
	// Here we only need deal with the klystron states.
	dfs_trace(dfs_trace::level_2)<<"State "<<nes<<": adjusting klystrons by "<<es.dEkly<<endl;
	if(es.dEkly!=0) {
		for(size_t k=0; k<theKlystrons.size(); k++) {
			Complex v = defkvals[k]*(1.0+es.dEkly);
			theKlystrons[k]->SetVoltagePhasor(v);
		}
	}
	else
		RestoreKlystrons(); // nominal gradient state
}

ConstantGradientAdjustment::ConstantGradientAdjustment()
: EnergyAdjustmentPolicy(), estates()
{
	estates.reserve(4);
	// Nominal state is always state 0
	estates.push_back(EnergyState(0,0));
}

ConstantGradientAdjustment::~ConstantGradientAdjustment()
{
	// nothing to do
}


