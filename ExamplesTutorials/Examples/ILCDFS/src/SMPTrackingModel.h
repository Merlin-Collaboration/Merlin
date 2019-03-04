/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_SMPTrackingModel
#define _h_SMPTrackingModel 1

#include "BeamDynamicsModel.h"
//#include "SMPBunch.h"
#include "SMPTracker.h"

class AcceleratorModel;
class BeamData;

namespace SMPTracking
{
class SMPBunch;
class WakeFieldProcess;
}

class SMPTrackingModel: public BeamDynamicsModel
{
public:

	// Constructor taking the number of slices (ns) and the
	// numner of macro-particles per slice (nps) to be used for
	// the bunch representation.
	explicit SMPTrackingModel(size_t ns, size_t nps);
	~SMPTrackingModel();

	void SetBeamline(const AcceleratorModel::Beamline& bl);
	void SetInitialBunch(const Bunch* bunch0);
	SMPTracking::SMPBunch* TrackBunch();
	void TrackThisBunch(Bunch* b);
	SMPTracking::SMPBunch* CreateBunch(const BeamData& beam0);

	void IncludeTransverseWakefield(bool flg);

private:

	size_t ns;
	size_t nps;
	SMPTracking::SMPTracker* tracker;
	SMPTracking::WakeFieldProcess* wfp;
	const SMPTracking::SMPBunch* cBunch0;
};

#endif
