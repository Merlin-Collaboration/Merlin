/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:19:43 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/StdField/TransRFfield.h"

Vector3D TransverseRFfield::GetBFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,0);
}

Vector3D TransverseRFfield::GetEFieldAt (const Point3D& x, double t) const
{
    double Er = TransverseRFfield::Er(x.z,t);
    return Vector3D(Er*cos(theta),Er*sin(theta),0);
}

Vector3D TransverseRFfield::GetForceAt (const Point3D& x, const Vector3D& v, double q, double t) const
{
   double Er = TransverseRFfield::Er(x.z,t);
   double Ex = Er*cos(theta);
   double Ey = Er*sin(theta);
   return Vector3D(q*Ex,q*Ey,0);
}

