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

// SWRFfield
#include "AcceleratorModel/StdField/SWRFfield.h"

Vector3D SWRFfield::GetBFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,0);
}

Vector3D SWRFfield::GetEFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,SWRFfield::Ez(x.z,t));
}

Vector3D SWRFfield::GetForceAt (const Point3D& x, const Vector3D& v, double q, double t) const
{
    return Vector3D(0,0,q*SWRFfield::Ez(x.z,t));
}

