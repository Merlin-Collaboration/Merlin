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

#ifndef GeneralRotation_h
#define GeneralRotation_h 1

#include "merlin_config.h"
// LinearAlgebra
#include "TLAS/LinearAlgebra.h"
// Space3D
#include "EuclideanGeometry/Space3D.h"
// RotationType
#include "EuclideanGeometry/RotationType.h"
// Rot3Drep
#include "EuclideanGeometry/Rot3Drep.h"
// AxisRotations
#include "EuclideanGeometry/AxisRotations.h"

//	A general (arbirary) 3-D rotation. Note that a general
//	rotation can represent a pure axis rotation.

class GeneralRotation : public Rot3Drep
{
public:

    //	Copy constructor.
    GeneralRotation (const GeneralRotation& gr);

    //	Contruction from a rotation matrix
    explicit GeneralRotation (const RealMatrix& r);

    //	Return an inverted rotation.
    virtual Rot3Drep* inv () const;

    //	Rotate the specified point and return the result.
    virtual Point3D rotate (const Point3D& p) const;

    //	Rotate the specified vector and return the result.
    virtual Vector3D rotate (const Vector3D& v) const;

    //	Dot this rotation with r.
    virtual Rot3Drep* dot (const Rot3Drep& r) const;

    //	Dot this rotation by a pure x rotation.
    virtual Rot3Drep* dotBy (const RotationX& rx) const;

    //	Dot this rotation by a pure y rotation.
    virtual Rot3Drep* dotBy (const RotationY& ry) const;

    //	Dot this rotation by a pure z rotation.
    virtual Rot3Drep* dotBy (const RotationZ& rz) const;

    //	Dot this rotation by a general rotation.
    virtual Rot3Drep* dotBy (const GeneralRotation& r) const;

    //	Return the type of rotation. Can be identity, xrot,
    //	yrot, zrot or general.
    virtual RotationType type () const;

    //	Returns true if this is a null rotation.
    virtual bool isIdentity () const;

    //	Returns true if a pure rotation about the x-axis.
    virtual bool isXrot () const;

    //	Returns true if a pure rotation about the y-axis.
    virtual bool isYrot () const;

    //	Returns true if a pure rotation about the z-axis.
    virtual bool isZrot () const;

    //	Rotation by 180 degrees about the X-axis.
    virtual Rot3Drep* rotXbyPI () const;

    //	Rotation by 180 degrees about the Y-axis.
    virtual Rot3Drep* rotYbyPI () const;

    //	Rotation by 180 degrees about the Z-axis.
    virtual Rot3Drep* rotZbyPI () const;

    //	Returns in m the 3x3 rotation matrix.
    virtual RealMatrix& getMatrix (RealMatrix& m) const;

private:

    //	Special constructors from two succesive axis rotations.
    GeneralRotation (const RotationX& rx, const RotationY& ry);
    GeneralRotation (const RotationX& rx, const RotationZ& rz);
    GeneralRotation (const RotationY& ry, const RotationX& rx);
    GeneralRotation (const RotationY& ry, const RotationZ& rz);
    GeneralRotation (const RotationZ& rz, const RotationX& rx);
    GeneralRotation (const RotationZ& rz, const RotationY& ry);
 
    //	Special constructors formed by the dot product of a
    //	general rotation with a pure axis rotation.
    GeneralRotation (const GeneralRotation& r, const RotationX& rx);
    GeneralRotation (const GeneralRotation& r, const RotationY& ry);
    GeneralRotation (const GeneralRotation& r, const RotationZ& rz);

    //	Default constructor.
    GeneralRotation ();

    RealMatrix m;

    friend class RotationX;
    friend class RotationY;
    friend class RotationZ;
};

inline GeneralRotation::GeneralRotation ()
        : m(IdentityMatrix(3))
{}

inline GeneralRotation::GeneralRotation (const RealMatrix& r)
        : m(r)
{}

#endif
