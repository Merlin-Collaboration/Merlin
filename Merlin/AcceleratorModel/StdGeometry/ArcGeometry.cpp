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

#include "AcceleratorModel/StdGeometry/ArcGeometry.h"

namespace {

inline Transform3D MakeTransform(double phi, double h)
{
    return Transform3D(Point3D((cos(phi)-1)/h,0,sin(phi)/h),Rotation3D::rotationY(-phi));
}

}; // end annonymous namespace


Transform3D ArcGeometry::GetGeometryTransform (double s0, double s) const throw (BeyondExtent)
{
    CheckBounds(s,s0);
    return MakeTransform((s-s0)*h,h);
}

Transform3D ArcGeometry::GetGeometryTransform (BoundaryPlane p) const
{
    double phi = p==entrance ? -GetAngle()/2 : GetAngle()/2;
    return MakeTransform(phi,h);
}

Transform3D ArcGeometry::GetTotalGeometryTransform () const
{
    return MakeTransform(GetAngle(),h);
}

