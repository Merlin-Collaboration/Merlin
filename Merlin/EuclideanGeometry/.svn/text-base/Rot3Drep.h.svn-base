/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef Rot3Drep_h
#define Rot3Drep_h 1

#include "merlin_config.h"
#include "TLAS/LinearAlgebra.h"
// Space3D
#include "EuclideanGeometry/Space3D.h"
// RotationType
#include "EuclideanGeometry/RotationType.h"

class IdentityRotation;
class GeneralRotation;
class RotationZ;
class RotationY;
class RotationX;

//	The abstract representation (letter class) for
//	Rotation3D.

class Rot3Drep
{
public:

    //	Virtual destructor.
    virtual ~Rot3Drep ();

    //	Rotate the specified point and return the result.
    virtual Point3D rotate (const Point3D& x) const = 0;

    //	Rotate the specified vector and return the result.
    virtual Vector3D rotate (const Vector3D& v) const = 0;

    //	Dot this rotation with r.
    virtual Rot3Drep* dot (const Rot3Drep& r) const = 0;

    //	Dot this rotation by a pure x rotation.
    virtual Rot3Drep* dotBy (const RotationX& rx) const = 0;

    //	Dot this rotation by a pure y rotation.
    virtual Rot3Drep* dotBy (const RotationY& ry) const = 0;

    //	Dot this rotation by a pure z rotation.
    virtual Rot3Drep* dotBy (const RotationZ& rz) const = 0;

    //	Dot this rotation by a general rotation.
    virtual Rot3Drep* dotBy (const GeneralRotation& r) const = 0;

    //	Return an inverted rotation.
    virtual Rot3Drep* inv () const = 0;

    //	Return the type of rotation. Can be identity, xrot,
    //	yrot, zrot or general.
    virtual RotationType type () const = 0;

    //	Returns true if this is a null rotation.
    virtual bool isIdentity () const;

    //	Returns true if a pure rotation about the x-axis.
    virtual bool isXrot () const;

    //	Returns true if a pure rotation about the y-axis.
    virtual bool isYrot () const;

    //	Returns true if a pure rotation about the z-axis.
    virtual bool isZrot () const;

    //	Static factory method. Constructs an identity (null)
    //	rotation.
    static Rot3Drep* identity ();

    //	Static factory method. Constructs a pure x rotation of
    //	angle radians.
    static Rot3Drep* rotationX (double angle);

    //	Static factory method. Constructs a pure y rotation of
    //	angle radians.
    static Rot3Drep* rotationY (double angle);

    //	Static factory method. Constructs a pure z rotation of
    //	angle radians.
    static Rot3Drep* rotationZ (double angle);

    //	Rotation by 180 degrees about the X-axis.
    virtual Rot3Drep* rotXbyPI () const = 0;

    //	Rotation by 180 degrees about the Y-axis.
    virtual Rot3Drep* rotYbyPI () const = 0;

    //	Rotation by 180 degrees about the Z-axis.
    virtual Rot3Drep* rotZbyPI () const = 0;

    //	Returns in m the 3x3 rotation matrix.
    virtual RealMatrix& getMatrix (RealMatrix& m) const = 0;

protected:

    //	Protected default constructor.
    Rot3Drep ();

private:

    //	Reference count for memory management.
    int refc;

    friend class Rotation3D;
};

#endif
