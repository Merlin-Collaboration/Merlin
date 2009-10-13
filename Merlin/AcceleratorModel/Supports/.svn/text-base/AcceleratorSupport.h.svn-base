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

#ifndef AcceleratorSupport_h
#define AcceleratorSupport_h 1

#include "merlin_config.h"
#include <vector>
#include "EuclideanGeometry/Space2D.h"
#include "EuclideanGeometry/Space3D.h"

//	Represents a single support structure. Support can be
//	translated in either x, y or z directions. Support is
//	primarily intended for ground motion application.

//	An array of AcceleratorSupport pointers.
class AcceleratorSupport;
typedef std::vector<AcceleratorSupport*> AcceleratorSupportList;

class AcceleratorSupport
{
public:

    AcceleratorSupport ();

    //	Sets the location of the support in the global
    //	coordinate frame. (x,z) represent the location in the
    //	accelerator plane, while s is the arc position of the
    //	support.
    void SetPosition (double s, double x, double z);

    //	Returns the arc position.
    double GetArcPosition () const;

    //	Returns the location of the support in the accelerator
    //	plane (x,z). Note that  Point2D::y here refers to the
    //	z-coordinate.
    Point2D GetLocation () const;

    //	Returns the linear distance from this support to a
    //	Support.
    double DistanceTo (const AcceleratorSupport& aSupport) const;

    //	Set the support offset to (x,y,z).
    void SetOffset (double x, double y, double z);

    //	Set the support offset to X
    void SetOffset (const Vector3D& X);

    //	Returns the current offset.
    const Vector3D& GetOffset () const;

    //	Increment the current offset by (dx,dy,dz),
    const Vector3D& IncrementOffset (double dx, double dy, double dz);

    //	Increment the current offset by dX.
    const Vector3D& IncrementOffset (const Vector3D& dX);

    //	Reset the offset to (0,0,0).
    void Reset ();

private:

    Vector3D offset;
    Point2D pos;
    double s_pos;
    bool modified;

    friend class SupportStructure;
};

inline AcceleratorSupport::AcceleratorSupport ()
        :offset(0,0,0),modified(false),pos(0,0),s_pos(0)
{}

inline void AcceleratorSupport::SetPosition (double s, double x, double z)
{
    pos.x=x;
    pos.y=z;
    s_pos=s;
}

inline double AcceleratorSupport::GetArcPosition () const
{
    return s_pos;
}

inline Point2D AcceleratorSupport::GetLocation () const
{
    return pos;
}

inline void AcceleratorSupport::SetOffset (double x, double y, double z)
{
    SetOffset(Vector3D(x,y,z));
}

inline void AcceleratorSupport::SetOffset (const Vector3D& X)
{
    offset=X;
    modified=true;
}

inline const Vector3D& AcceleratorSupport::GetOffset () const
{
    return offset;
}

inline const Vector3D& AcceleratorSupport::IncrementOffset (double dx, double dy, double dz)
{
    return IncrementOffset(Vector3D(dx,dy,dz));
}

inline const Vector3D& AcceleratorSupport::IncrementOffset (const Vector3D& dX)
{
    modified=true;
    return offset+=dX;
}

inline void AcceleratorSupport::Reset ()
{
    modified=true;
    offset=Vector3D(0,0,0);
}

#endif
