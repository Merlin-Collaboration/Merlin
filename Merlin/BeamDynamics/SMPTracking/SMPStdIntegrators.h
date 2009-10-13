/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (1999)
* 
* file Merlin\BeamDynamics\SliceMPTracking\StdIntegrators.h
* last modified 11/12/01 15:32:02
*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
* Copyright (c) 2001 by The Merlin Collaboration.
* - ALL RIGHTS RESERVED - 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all coIes and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/

#ifndef SliceMPTracking_StdIntegrators_h
#define SliceMPTracking_StdIntegrators_h 1


#include "merlin_config.h"
#include "BeamDynamics/SMPTracking/SMPComponentTracker.h"


// SectorBend
#include "AcceleratorModel/StdComponent/SectorBend.h"
// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"
// SWRFStructure
#include "AcceleratorModel/StdComponent/SWRFStructure.h"
// TWRFStructure
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
// Drift
#include "AcceleratorModel/StdComponent/Drift.h"
// Marker
#include "AcceleratorModel/StdComponent/Marker.h"
// MatrixMaps
#include "BasicTransport/MatrixMaps.h"
// TransportMatrix
#include "BasicTransport/TransportMatrix.h"
// Monitor
#include "AcceleratorModel/StdComponent/Monitor.h"
// Solenoid
#include "AcceleratorModel/StdComponent/Solenoid.h"


namespace SMPTracking {

class DriftCI : public SMPComponentTracker::Integrator<Drift> {
protected:
    void TrackStep (double ds);
};

class TWRFStructureCI : public SMPComponentTracker::Integrator<TWRFStructure> {
protected:
    void TrackStep (double ds);
    void TrackEntrance();
    void TrackExit();
private:
    void ApplyEndField(double);
};

class SectorBendCI : public SMPComponentTracker::Integrator<SectorBend> {
protected:
    void TrackStep (double ds);
    void TrackEntrance();
    void TrackExit();

    //	Used to apply a linear pole face rotation to the current
    //	bunch.
    void ApplyPoleFaceRotation (double h, const SectorBend::PoleFace& pf);
};

class RectMultipoleCI : public SMPComponentTracker::Integrator<RectMultipole> {
protected:
    void TrackStep (double ds);
};

class MonitorCI : public SMPComponentTracker::Integrator<Monitor>  {
protected:
    void TrackStep (double ds);
};

class SWRFStructureCI : public SMPComponentTracker::Integrator<SWRFStructure> {
protected:
    void TrackStep (double ds);
};

class SolenoidCI : public SMPComponentTracker::Integrator<Solenoid> {
protected:
    void TrackStep (double ds);
};

class MarkerCI : public SMPComponentTracker::Integrator<Marker> {
protected:
    void TrackStep (double ds) { return; }
};
};

#endif
