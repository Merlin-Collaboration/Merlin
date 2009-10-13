/////////////////////////////////////////////////////////////////////////
// class ParticleBunchModel
// Uses a particle ensemble to represent the bunch. 
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

#ifndef _h_ParticleTrackingModel
#define _h_ParticleTrackingModel 1

#include "BeamDynamicsModel.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"

class AcceleratorModel;
class BeamData;

class ParticleTrackingModel : public BeamDynamicsModel {	
public:

	explicit ParticleTrackingModel(size_t npart);
	~ParticleTrackingModel();

	void SetBeamline(const AcceleratorModel::Beamline& bl);	
	void SetInitialBunch(const Bunch* bunch0);	
	ParticleTracking::ParticleBunch* TrackBunch();	
	void TrackThisBunch(Bunch* b);
	ParticleTracking::ParticleBunch* CreateBunch(const BeamData& beam0);

	void IncludeTransverseWakefield(bool flg);	

private:

	size_t np;
	ParticleTracking::ParticleTracker* tracker;
	const ParticleTracking::ParticleBunch* cBunch0;
};

#endif
