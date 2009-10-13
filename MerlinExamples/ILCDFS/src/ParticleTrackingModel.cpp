#include "ParticleTrackingModel.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "ILCDFS_IO.h"

using namespace ParticleTracking;

ParticleTrackingModel::ParticleTrackingModel(size_t npart) 
: BeamDynamicsModel("PARTICLE TRACKING"),np(npart),cBunch0(0)
{
	dfs_trace(dfs_trace::level_1)<<"PARTICLE TRACKING initialised with "<<npart<<" particles"<<endl;
	tracker = new ParticleTracker();
}

ParticleTrackingModel::~ParticleTrackingModel()
{
	if(tracker)
		delete tracker;
}

void ParticleTrackingModel::SetBeamline(const AcceleratorModel::Beamline& bl)
{
	tracker->SetOutput(output);
	tracker->SetBeamline(bl);
}

void ParticleTrackingModel::SetInitialBunch(const Bunch* bunch0)
{
	// We assume here that the bunch is a ParticleBunch
	cBunch0 = static_cast<const ParticleBunch*>(bunch0);
}

ParticleBunch* ParticleTrackingModel::TrackBunch()
{
	ParticleBunch* b = new ParticleBunch(*cBunch0);
	tracker->SetOutput(output);
	tracker->Track(b);
	return b;
}

void ParticleTrackingModel::TrackThisBunch(Bunch* b)
{
	tracker->Track(static_cast<ParticleBunch*>(b));
}

ParticleBunch* ParticleTrackingModel::CreateBunch(const BeamData& beam0)
{
	return ParticleBunchConstructor(beam0,np).ConstructParticleBunch();
}

void ParticleTrackingModel::IncludeTransverseWakefield(bool flg)
{
	// TODO Auto-generated method stub
}