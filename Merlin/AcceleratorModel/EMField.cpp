/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/EMField.h"

EMField::~EMField ()
{
    // Nothing to do
}

Vector3D EMField::GetForceAt (const Point3D& x, const Vector3D& v, double q, double t) const
{
    Vector3D B = GetBFieldAt(x,t);
    Vector3D E = GetEFieldAt(x,t);
    return q*(E+cross(v,B));
}

