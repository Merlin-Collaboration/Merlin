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

// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int RectMultipole::ID = UniqueIndex();

RectMultipole::RectMultipole (const string& id, double length, int npole, double b, double r0, bool skew)
        : TAccCompGF<RectangularGeometry,MultipoleField>(id,new RectangularGeometry(length),
                new MultipoleField(npole,b,r0,skew))
{}

RectMultipole::RectMultipole (const string& id, double len, int np, double b, bool skew)
        : TAccCompGF<RectangularGeometry,MultipoleField>(id,new RectangularGeometry(len),
                new MultipoleField(np,b,skew))
{}

int RectMultipole::GetIndex () const
{
    return ID;
}

void RectMultipole::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

void RectMultipole::RotateY180 ()
{
    GetField().RotateY180();
}

