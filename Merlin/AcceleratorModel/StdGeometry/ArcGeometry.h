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

#ifndef ArcGeometry_h
#define ArcGeometry_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdGeometry/CenteredGeometry.h"

//	An ArcGeometry represents a constant radius curve in the
//	local x-z plane. By convention, a positive curvature (h)
//	defines a curve towards negative x. Transformations
//	between two points (s1,s2) on an ArcGeometry are
//	specified by a translation in the x-z plane of
//	[(cos(phi)-1)/h, 0, sin(phi)/h], and a rotation about
//	the y-axis of -phi (phi = h*(s2-s1)).
//
//	An ArcGeometry can have an additional tilt, which is
//	defined as a rotation about the local z-axis of the
//	entrance plane by an angle theta, with an additional
//	rotation about the exit plane z-axis by -theta.

class ArcGeometry : public CenteredGeometry
{
public:

    //	Constructor taking the arc length and the curvature
    //	(1/r) of the geometry.
    ArcGeometry (double l, double curv);

    //	Return the curvature of the ArcGeometry.
    double GetCurvature () const;

    //	Return the total arc angle of the geometry.
    double GetAngle () const;

    //	Returns the arc transform from s0 to s (angle=h*(s-s0)).
    virtual Transform3D GetGeometryTransform (double s0, double s) const throw (BeyondExtent);

    //	Returns the arc transform for either -angle/2 or
    //	+angle/2 for the entrance and exit planes respectively.
    virtual Transform3D GetGeometryTransform (BoundaryPlane p) const;

    //	Returns the arc transformation from the entrance plane
    //	frame to the exit plane frame.
    virtual Transform3D GetTotalGeometryTransform () const;

    //	Sets the curvature (=1/r) of the arc.
    void SetCurvature (double curv);

    //	Sets the tilt of the geomety in radians.
    void SetTilt (double t);

    //	Returns the tilt of the geometry in radians.
    double GetTilt () const;

private:

    double h;
    double tilt;
};

inline ArcGeometry::ArcGeometry (double l, double curv)
        : CenteredGeometry(l),h(curv),tilt(0)
{}

inline double ArcGeometry::GetCurvature () const
{
    return h;
}

inline double ArcGeometry::GetAngle () const
{
    return len*h;
}

inline void ArcGeometry::SetCurvature (double curv)
{
    h=curv;
}

inline void ArcGeometry::SetTilt (double t)
{
    tilt=t;
}

inline double ArcGeometry::GetTilt () const
{
    return tilt;
}

#endif
