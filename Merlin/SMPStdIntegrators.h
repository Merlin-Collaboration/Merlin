/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SliceMPTracking_StdIntegrators_h
#define SliceMPTracking_StdIntegrators_h 1

#include "merlin_config.h"
#include "SMPComponentTracker.h"

#include "SectorBend.h"
#include "RectMultipole.h"
#include "SWRFStructure.h"
#include "TWRFStructure.h"
#include "Drift.h"
#include "Marker.h"
#include "MatrixMaps.h"
#include "TransportMatrix.h"
#include "Monitor.h"
#include "Solenoid.h"

namespace SMPTracking
{

class DriftCI: public SMPComponentTracker::Integrator<Drift>
{
protected:
	void TrackStep(double ds);
};

class TWRFStructureCI: public SMPComponentTracker::Integrator<TWRFStructure>
{
protected:
	void TrackStep(double ds);
	void TrackEntrance();
	void TrackExit();
private:
	void ApplyEndField(double);
};

class SectorBendCI: public SMPComponentTracker::Integrator<SectorBend>
{
protected:
	void TrackStep(double ds);
	void TrackEntrance();
	void TrackExit();

	/**
	 *	Used to apply a linear pole face rotation to the current
	 *	bunch.
	 */
	void ApplyPoleFaceRotation(double h, const SectorBend::PoleFace& pf);
};

class RectMultipoleCI: public SMPComponentTracker::Integrator<RectMultipole>
{
protected:
	void TrackStep(double ds);
};

class MonitorCI: public SMPComponentTracker::Integrator<Monitor>
{
protected:
	void TrackStep(double ds);
};

class SWRFStructureCI: public SMPComponentTracker::Integrator<SWRFStructure>
{
protected:
	void TrackStep(double ds);
};

class SolenoidCI: public SMPComponentTracker::Integrator<Solenoid>
{
protected:
	void TrackStep(double ds);
};

class MarkerCI: public SMPComponentTracker::Integrator<Marker>
{
protected:
	void TrackStep(double ds)
	{
		return;
	}
};
}

#endif
