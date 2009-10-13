/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/26 20:12:15 $
// $Revision: 1.6 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ParticleTracking_StdIntegrators_h
#define ParticleTracking_StdIntegrators_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdComponent/SectorBend.h"
#include "AcceleratorModel/StdComponent/RectMultipole.h"
#include "AcceleratorModel/StdComponent/SWRFStructure.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "AcceleratorModel/StdComponent/Drift.h"
#include "AcceleratorModel/StdComponent/Marker.h"
#include "AcceleratorModel/StdComponent/Monitor.h"
#include "AcceleratorModel/StdComponent/Solenoid.h"

#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"
#include "BeamDynamics/ParticleTracking/Integrators/ParticleMapPI.h"

#define DECL_SIMPLE_INTG(I,C) class I : \
	public ParticleComponentTracker::Integrator< C > { \
	public: void TrackStep(double); };

namespace ParticleTracking {

// common integrators
DECL_SIMPLE_INTG(MonitorCI,Monitor);
DECL_SIMPLE_INTG(MarkerCI,Marker);
DECL_SIMPLE_INTG(SolenoidCI,Solenoid);
//DECL_SIMPLE_INTG(ParticleMapCI,ParticleMapComponent);

namespace THIN_LENS {

DECL_SIMPLE_INTG(DriftCI,Drift);
DECL_SIMPLE_INTG(RectMultipoleCI,RectMultipole);
DECL_SIMPLE_INTG(TWRFStructureCI,TWRFStructure);
DECL_SIMPLE_INTG(SWRFStructureCI,SWRFStructure);

class SectorBendCI : public ParticleComponentTracker::Integrator<SectorBend> {
public:
    void TrackStep(double);
    void TrackEntrance();
    void TrackExit();
protected:
    void ApplyPoleFaceRotation(double h, const SectorBend::PoleFace& pf);
};

DECL_INTG_SET(ParticleComponentTracker,StdISet)
};

namespace TRANSPORT {

DECL_SIMPLE_INTG(DriftCI,Drift);
DECL_SIMPLE_INTG(RectMultipoleCI,RectMultipole);

class SectorBendCI : public ParticleComponentTracker::Integrator<SectorBend> {
public:
    void TrackStep(double);
    void TrackEntrance();
    void TrackExit();
protected:
    void ApplyPoleFaceRotation(const SectorBend::PoleFace* pf);
};

DECL_INTG_SET(ParticleComponentTracker,StdISet)
};

}; // end namespace ParticleTracking

#endif
