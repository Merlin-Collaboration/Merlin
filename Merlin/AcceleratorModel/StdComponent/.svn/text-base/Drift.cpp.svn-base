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

// Drift
#include "AcceleratorModel/StdComponent/Drift.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int Drift::ID = UniqueIndex();

Drift::Drift (const string& id, double len)
        : TAccCompG<RectangularGeometry>(id,new RectangularGeometry(len))
{}

const string& Drift::GetType () const
{
    _TYPESTR(Drift);
}

int Drift::GetIndex () const
{
    return  ID;
}

void Drift::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

void Drift::RotateY180 ()
{
    // nothing to do
}

ModelElement* Drift::Copy () const
{
    return new Drift(*this);
}

