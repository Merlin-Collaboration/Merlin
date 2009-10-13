/////////////////////////////////////////////////////////////////////////
// class SMPTrackingModel
// Uses Sliced-Macro-Particles to represent a bunch 
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

#ifndef _h_SMPTrackingModel
#define _h_SMPTrackingModel 1

#include "BeamDynamicsModel.h"
//#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamDynamics/SMPTracking/SMPTracker.h"

class AcceleratorModel;
class BeamData;

namespace SMPTracking {
	class SMPBunch;
	class WakeFieldProcess;
};

class SMPTrackingModel : public BeamDynamicsModel {	
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
