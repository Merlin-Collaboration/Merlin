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

#ifndef AcceleratorGeometry_h
#define AcceleratorGeometry_h 1

#include "merlin_config.h"
#include <utility>

// Transform3D
#include "EuclideanGeometry/Transform3D.h"

//	Represents a frame of reference for a section of
//	accelerator lattice. An AcceleratorGeometry can be
//	considered a type of three-dimensional space line
//	(R(s)), which is characterised by a single scalar s, the
//	distance along the space line from the origin. Each
//	AcceleratorGeometry has a specific length which bounds
//	the allowed values of s (with respect to the local
//	geometry origin).  At each  position s on the geometry,
//	a local rectangular coordinate frame can be uniquely
//	defined, with its origin at the point s, and its z-axis
//	tangential to the geometry at s. The orientation of the
//	local x- and y-axes are also uniquely determined by the
//	sum of any rotations applied going from the origin to s.
//
//	The primary responsibilty for an AcceleratorGeometry
//	object is to supply transformations between coordinate
//	frames defined on that geometry.


class AcceleratorGeometry
{
public:
    //	Exception thrown indicating that an s-distance was
    //	outside of the current geometry extent.
    class BeyondExtent {};

    typedef std::pair<double,double> Extent;

    //	A BoundaryPlane is the X-Y plane (z=0) of the coordinate
    //	frame defined at the entrance (start) or exit (end) of
    //	the Geometry.
    typedef enum {entrance,exit} BoundaryPlane;

public:
    //	Virtual destructor.
    virtual ~AcceleratorGeometry ();

    //	Return the three-dimensional transformation from the
    //	frame at s0 to the frame at s. s and s0 are in the
    //	geometry's s-frame, and must be within the geometry
    //	extents.
    virtual Transform3D GetGeometryTransform (double s0, double s) const throw (BeyondExtent) = 0;

    //	Return the three-dimensional transformation from the
    //	local origin to the frame at s. s is  in the geometry's
    //	s-frame, and must be within the geometry extents.
    Transform3D GetGeometryTransform (double s) const throw (BeyondExtent);

    //	Returns the transformation from the geometry origin to
    //	the specified boundary plane.
    virtual Transform3D GetGeometryTransform (BoundaryPlane p) const = 0;

    //	Returns the transformation from the entrance plane frame
    //	to the exit plane frame.
    virtual Transform3D GetTotalGeometryTransform () const;

    //	Returns the local extent of this geometry.
    virtual Extent GetGeometryExtent () const = 0;

    //	Returns the total arc-length of the geometry.
    virtual double GetGeometryLength () const = 0;
};

inline AcceleratorGeometry::~AcceleratorGeometry ()
{}

inline Transform3D AcceleratorGeometry::GetGeometryTransform (double s) const throw (BeyondExtent)
{
    return GetGeometryTransform(0,s);
}

inline double ToLength(const AcceleratorGeometry::Extent& extent)
{
    return extent.second-extent.first;
}


#endif
