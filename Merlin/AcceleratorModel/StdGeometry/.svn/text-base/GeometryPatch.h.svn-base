/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:17:37 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef GeometryPatch_h
#define GeometryPatch_h 1

#include "merlin_config.h"
#include <utility>

// Transform3D
#include "EuclideanGeometry/Transformable.h"
#include "AcceleratorModel/AcceleratorGeometry.h"

// A geometry patch, representing an arbitrary
// transformation of the accelerator geometry.
// GeometryPatch objects have zero extents.
//
// Note: this class was developed primarilly to
//       support MAD-like SROT elements.

class GeometryPatch : public AcceleratorGeometry, public Transformable
{
public:

    virtual Transform3D GetGeometryTransform (double s0, double s) const throw (BeyondExtent)
    {
        if(s0!=0 || s!=0)
            throw BeyondExtent();
        return local_T ? *local_T : Transform3D();
    }

    virtual Transform3D GetGeometryTransform (BoundaryPlane p) const
    {
        return (p==entrance && local_T) ? *local_T : Transform3D();
    }

    virtual Transform3D GetTotalGeometryTransform () const
    {
        return local_T ? *local_T : Transform3D();
    }

    //	Returns the local extent of this geometry.
    virtual Extent GetGeometryExtent () const { return Extent(0,0); }

    //	Returns the total arc-length of the geometry.
    virtual double GetGeometryLength () const { return 0.0; }
};

#endif
