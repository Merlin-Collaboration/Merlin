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

// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// Marker
#include "AcceleratorModel/StdComponent/Marker.h"

const int Marker::ID = UniqueIndex();

void Marker::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

int Marker::GetIndex () const
{
    return ID;
}

const string& Marker::GetType () const
{
    _TYPESTR(Marker)
}

ModelElement* Marker::Copy () const
{
    return new Marker(*this);
}

void Marker::RotateY180 ()
{
    // Nothing to do
}

