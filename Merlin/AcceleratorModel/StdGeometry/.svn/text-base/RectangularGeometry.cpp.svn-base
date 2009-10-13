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

#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"

Transform3D RectangularGeometry::GetGeometryTransform (double s0, double s) const throw (BeyondExtent)
{
    CheckBounds(s,s0);
    return Transform3D::translation(0,0,s-s0);
}

Transform3D RectangularGeometry::GetGeometryTransform (BoundaryPlane p) const
{
    double s = p==entrance ? -len/2 : len/2;
    return Transform3D::translation(0,0,s);
}

Transform3D RectangularGeometry::GetTotalGeometryTransform () const
{
    return Transform3D::translation(0,0,len);
}

