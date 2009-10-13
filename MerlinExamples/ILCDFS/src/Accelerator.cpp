/////////////////////////////////////////////////////////////////////////
// Class Accelerator implementation
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

// Merlinlib 
#include "Channels/Channels.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/Bunch.h"
#include "BeamModel/BeamData.h"

// ILCDFS
#include "Accelerator.h"
#include "BeamDynamicsModel.h"
#include "ILCDFS_IO.h"

namespace {

	// Function used for sorting the klystron array
	bool KlystronSort(const Klystron* k1, const Klystron* k2)
	{
		vector<size_t> i;
		k1->GetBeamlineIndecies(i);
		size_t i1 = i.front();
		i.clear();
		k2->GetBeamlineIndecies(i);
		size_t i2 = i.front();
		return i1<i2;
	}

}; // end of anonymous namespace

Accelerator::Accelerator(const std::string& name, AcceleratorModel* aModel, BeamData* ibeamdat)
: itsName(name),itsAccModel(aModel),itsTracker(0),beam0(ibeamdat),cachedBunches(),
currentSegment(0,0),allowIncrTracking(false)
{}

Accelerator::~Accelerator()
{
	if(itsAccModel)
		delete itsAccModel;
	if(itsTracker)
		delete itsTracker;
	if(beam0)
		delete beam0;
}

const std::string& Accelerator::GetName() const
{
	return itsName;
}

void Accelerator::SetBeamDynamicsModel(BeamDynamicsModel* bdm)
{
	itsTracker = bdm;
	cachedBunches.clear();
	dfs_trace(dfs_trace::level_1)<<itsName<<" using "<<itsTracker->GetName()<<endl;
}

void Accelerator::AllowIncrementalTracking(bool flg)
{
	dfs_trace(dfs_trace::level_1)<<itsName<<" using incremental tracking = ";
	if(flg)
		dfs_trace(dfs_trace::level_1)<<"YES"<<endl;
	else
		dfs_trace(dfs_trace::level_1)<<"NO"<<endl;


	allowIncrTracking = flg;
}

void Accelerator::SetActiveBeamlineSegment(const DFS_Segment& seg)
{
	dfs_trace(dfs_trace::level_3)<<itsName<<" active segment set to "<<seg<<endl;
	currentSegment = seg;
}

void Accelerator::TrackNewBunchThroughModel()
{
	AcceleratorModel::Beamline bline = itsAccModel->GetBeamline();
	itsTracker->SetBeamline(bline);
	Bunch* b = itsTracker->CreateBunch(*beam0);
	itsTracker->TrackThisBunch(b);
	dfs_trace(dfs_trace::level_3)<<"final energy = "<<b->GetReferenceMomentum()<<" GeV"<<endl;
	SetActiveBeamlineSegment(currentSegment);
	delete b;
}

void Accelerator::TrackBeam(size_t nstate)
{
	dfs_trace(dfs_trace::level_3)<<itsName<<" tracking bunch for state "<<nstate;

	CachedBunch& cb = cachedBunches[nstate];

	if(allowIncrTracking) {
		dfs_trace(dfs_trace::level_3)<<" using incremental tracking";
		// check to see if this bunch is already at the 
		// entrance to the current segment. Exception when we are at the beginning
		// of the beamline.
		if(currentSegment.first!=0 && cb.location+1!=currentSegment.first) {
			// Calculate indecies (not location==0 is a special case)
			size_t n1 = cb.location == 0 ? 0 : cb.location+1;
			size_t n2 = currentSegment.first-1;
			dfs_trace(dfs_trace::level_3)<<"\n  - incrementing beam "<<nstate<<" from "<<n1<<" to "<<n2;

			AcceleratorModel::Beamline bline = itsAccModel->GetBeamline(n1,n2);
			itsTracker->SetBeamline(bline);
			itsTracker->TrackThisBunch(cb.bunch); // update the bunch
			cb.location = n2;
		}
	}
	dfs_trace(dfs_trace::level_3)<<endl;

	// If we are not doing incremental tracking, we always
	// track from the beginning of the beamline
	size_t n1 = allowIncrTracking ? currentSegment.first : 0;
	size_t n2 = currentSegment.second;
	AcceleratorModel::Beamline bline = itsAccModel->GetBeamline(n1,n2);

	itsTracker->SetBeamline(bline);
	itsTracker->SetInitialBunch(cb.bunch);
	Bunch* rb = itsTracker->TrackBunch(); // tracks a copy of the current bunch.
	dfs_trace(dfs_trace::level_3)<<"final energy = "<<rb->GetReferenceMomentum()<<" GeV"<<endl;
	delete rb;
}

size_t Accelerator::GetMonitorChannels(Plane p, ROChannelArray& bpmChannels)
{
	AcceleratorModel::Beamline bline = 
		itsAccModel->GetBeamline(currentSegment.first,currentSegment.second);

	std::vector<ROChannel*> bc;

	if(p==x_only || p==x_and_y)
		itsAccModel->GetROChannels(bline,"BPM.*.X",bc);

	if(p==y_only || p==x_and_y)
		itsAccModel->GetROChannels(bline,"BPM.*.Y",bc);

	bpmChannels.SetChannels(bc);
	return bpmChannels.Size();
}

size_t Accelerator::GetCorrectorChannels(Plane p, RWChannelArray& corrChannels)
{
	AcceleratorModel::Beamline bline = 
		itsAccModel->GetBeamline(currentSegment.first,currentSegment.second);

	std::vector<RWChannel*> cc;

	if(p==x_only || p==x_and_y)
		itsAccModel->GetRWChannels(bline,"XCor.*.B0",cc);

	if(p==y_only || p==x_and_y)
		itsAccModel->GetRWChannels(bline,"YCor.*.B0",cc);

	corrChannels.SetChannels(cc);
	return corrChannels.Size();
}

size_t Accelerator::GetKlystrons(KlystronArray& klys)
{
	itsAccModel->ExtractTypedElements(klys);
	// There is no guarantee that the klystrons are
	// in beamline order, so we must sort them.
	std::sort(klys.begin(),klys.end(),KlystronSort);
	return klys.size();
}

void Accelerator::InitialiseTracking(size_t nstates, ReferenceParticleArray& refplist)
{
	cachedBunches.clear();
	refplist.clear();
	cachedBunches.reserve(nstates);
	refplist.reserve(nstates);

	while((nstates--)>0) {
		Bunch* b = itsTracker->CreateBunch(*beam0);
		refplist.push_back(b);
		cachedBunches.push_back(CachedBunch(b));
	}
}

size_t Accelerator::GetBeamlineIndecies(const std::string& cpat, IntegerArray& indecies) const
{
	indecies.clear();
	return itsAccModel->GetIndecies(cpat,indecies);
}

DFS_Segment Accelerator::GetBeamlineRange() const
{
	AcceleratorModel::Beamline bl = itsAccModel->GetBeamline();
	return DFS_Segment(bl.first_index(),bl.last_index());
}
