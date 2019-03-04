/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ParticleTrackingModel
#define _h_ParticleTrackingModel 1

#include "BeamDynamicsModel.h"
#include "ParticleBunch.h"
#include "ParticleTracker.h"

class AcceleratorModel;
class BeamData;

class ParticleTrackingModel: public BeamDynamicsModel
{
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
