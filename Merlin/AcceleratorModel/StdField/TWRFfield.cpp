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

#include "AcceleratorModel/StdField/TWRFfield.h"

Vector3D TWRFfield::GetBFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,0);
}

Vector3D TWRFfield::GetEFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,TWRFfield::Ez(x.z,t));
}

Vector3D TWRFfield::GetForceAt (const Point3D& x, const Vector3D& v, double q, double t) const
{
    return Vector3D(0,0,q*TWRFfield::Ez(x.z,t));
}

