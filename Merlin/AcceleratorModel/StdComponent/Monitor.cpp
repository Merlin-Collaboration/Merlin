/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

// Bunch
#include "BeamModel/Bunch.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// Monitor
#include "AcceleratorModel/StdComponent/Monitor.h"

bool Monitor::all_inactive = false;
const int Monitor::ID = UniqueIndex();

Monitor::Monitor (const string& id, double len, double mpt)
        : TAccCompG<RectangularGeometry>(id,new RectangularGeometry(len)),mp(mpt),active(true)
{}

Monitor::~Monitor ()
{}

void Monitor::MakeMeasurement (const Bunch& )
{}

void Monitor::SetMeasurementPt (double mpt) throw (AcceleratorGeometry::BeyondExtent)
{
    //GetGeometry().CheckBounds(mpt); // might throw
    mp=mpt;
}

double Monitor::GetMeasurementPt () const
{
    return mp;
}

void Monitor::RotateY180 ()
{
    reflected=!reflected;
}

int Monitor::GetIndex () const
{
    return ID;
}

const string& Monitor::GetType () const
{
    _TYPESTR(Monitor)
}

void Monitor::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent)
}

ModelElement* Monitor::Copy () const
{
    return new Monitor(*this);
}

