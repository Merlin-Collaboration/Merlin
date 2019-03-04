/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ParticleTracking_StdIntegrators_h
#define ParticleTracking_StdIntegrators_h 1

#include "merlin_config.h"
#include "SectorBend.h"
#include "RectMultipole.h"
#include "SWRFStructure.h"
#include "TWRFStructure.h"
#include "Drift.h"
#include "Marker.h"
#include "Monitor.h"
#include "Solenoid.h"

#include "ParticleComponentTracker.h"
#include "ParticleMapPI.h"

#define DECL_SIMPLE_INTG(I, C) class I: \
	public ParticleComponentTracker::Integrator<C> { \
	public: void TrackStep(double); };

namespace ParticleTracking
{

// common integrators
DECL_SIMPLE_INTG(MonitorCI, Monitor)
DECL_SIMPLE_INTG(MarkerCI, Marker)
DECL_SIMPLE_INTG(SolenoidCI, Solenoid)
//DECL_SIMPLE_INTG(ParticleMapCI,ParticleMapComponent)

namespace THIN_LENS
{

DECL_SIMPLE_INTG(DriftCI, Drift)
DECL_SIMPLE_INTG(RectMultipoleCI, RectMultipole)
DECL_SIMPLE_INTG(TWRFStructureCI, TWRFStructure)
DECL_SIMPLE_INTG(SWRFStructureCI, SWRFStructure)

class SectorBendCI: public ParticleComponentTracker::Integrator<SectorBend>
{
public:
	void TrackStep(double);
	void TrackEntrance();
	void TrackExit();
protected:
	void ApplyPoleFaceRotation(double h, const SectorBend::PoleFace& pf);
};

DECL_INTG_SET(ParticleComponentTracker, StdISet)
} //end namespace THIN_LENS

namespace TRANSPORT
{

DECL_SIMPLE_INTG(DriftCI, Drift)
DECL_SIMPLE_INTG(RectMultipoleCI, RectMultipole)

class SectorBendCI: public ParticleComponentTracker::Integrator<SectorBend>
{
public:
	void TrackStep(double);
	void TrackEntrance();
	void TrackExit();
protected:
	void ApplyPoleFaceRotation(const SectorBend::PoleFace* pf);
};

DECL_INTG_SET(ParticleComponentTracker, StdISet)
} //end namespace TRANSPORT

} // end namespace ParticleTracking

#endif
