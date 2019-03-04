/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ParticleTrackingModel.h"
#include "WakeFieldProcess.h"
#include "ILCDFS_IO.h"
#include "ParticleBunch.h"
#include "ParticleDistributionGenerator.h"

using namespace ParticleTracking;

ParticleTrackingModel::ParticleTrackingModel(size_t npart) :
	BeamDynamicsModel("PARTICLE TRACKING"), np(npart), cBunch0(nullptr)
{
	dfs_trace(dfs_trace::level_1) << "PARTICLE TRACKING initialised with " << npart << " particles" << endl;
	tracker = new ParticleTracker();
}

ParticleTrackingModel::~ParticleTrackingModel()
{
	if(tracker)
	{
		delete tracker;
	}
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
	return new ParticleBunch(np, NormalParticleDistributionGenerator(), beam0);
}

void ParticleTrackingModel::IncludeTransverseWakefield(bool flg)
{
	// TODO Auto-generated method stub
}
