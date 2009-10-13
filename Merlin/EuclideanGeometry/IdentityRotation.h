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

#ifndef IdentityRotation_h
#define IdentityRotation_h 1

#include "merlin_config.h"
// Space3D
#include "EuclideanGeometry/Space3D.h"
// RotationType
#include "EuclideanGeometry/RotationType.h"
// Rot3Drep
#include "EuclideanGeometry/Rot3Drep.h"

class GeneralRotation;
class RotationZ;
class RotationY;
class RotationX;

//	A null or identity rotation.

class IdentityRotation : public Rot3Drep
{
public:

    //	Return an inverted rotation.
    virtual Rot3Drep* inv () const;

    //	Rotate the specified point and return the result.
    virtual Point3D rotate (const Point3D& p) const;

    //	Rotate the specified vector and return the result.
    virtual Vector3D rotate (const Vector3D& v) const;

    //	Returns true.
    virtual bool isIdentity () const;

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

    //	Rotation by 180 degrees about the X-axis.
    virtual Rot3Drep* rotXbyPI () const;

    //	Rotation by 180 degrees about the Y-axis.
    virtual Rot3Drep* rotYbyPI () const;

    //	Rotation by 180 degrees about the Z-axis.
    virtual Rot3Drep* rotZbyPI () const;

    //	Returns in m the 3x3 rotation matrix.
    virtual RealMatrix& getMatrix (RealMatrix& m) const;
};

#endif
