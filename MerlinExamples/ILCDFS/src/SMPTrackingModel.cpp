#include "SMPTrackingModel.h"
#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamDynamics/SMPTracking/SMPBunchConstructor.h"
#include "BeamDynamics/SMPTracking/SMPWakeFieldProcess.h"
#include "ILCDFS_IO.h"

using namespace SMPTracking;

SMPTrackingModel::SMPTrackingModel(size_t nslice, size_t mpps) 
: BeamDynamicsModel("SMP TRACKING"),ns(nslice),nps(mpps),cBunch0(0)
{
	dfs_trace(dfs_trace::level_1)<<"SMP TRACKING initialised with n_slice = "<<ns<<", n_mp = "<<nps<<endl;
	tracker = new SMPTracker();
	wfp = new WakeFieldProcess(1);
	wfp->ApplyImpulseAt(WakeFieldProcess::atCentre);
	tracker->AddProcess(wfp);
}

SMPTrackingModel::~SMPTrackingModel()
{
	if(tracker)
		delete tracker;
}

void SMPTrackingModel::SetBeamline(const AcceleratorModel::Beamline& bl)
{
	tracker->SetBeamline(bl);
}

void SMPTrackingModel::SetInitialBunch(const Bunch* bunch0)
{
	// We assume here that the bunch is a SMPBunch
	cBunch0 = static_cast<const SMPBunch*>(bunch0);
}

SMPBunch* SMPTrackingModel::TrackBunch()
{
	SMPBunch* b = new SMPBunch(*cBunch0);
	tracker->SetOutput(output);
	tracker->Track(b);
	return b;
}

void SMPTrackingModel::TrackThisBunch(Bunch* b)
{
	tracker->SetOutput(output);
	tracker->Track(static_cast<SMPBunch*>(b));
}

SMPBunch* SMPTrackingModel::CreateBunch(const BeamData& beam0)
{
	return SMPBunchConstructor(beam0,ns,nps).ConstructSMPBunch();
}

void SMPTrackingModel::IncludeTransverseWakefield(bool flg)
{
	dfs_trace(dfs_trace::level_2)<<"including transverse wakes: "<<flg<<endl;
	wfp->IncludeTransverseWake(flg);
}
