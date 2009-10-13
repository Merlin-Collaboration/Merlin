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

#ifndef Rotation3D_h
#define Rotation3D_h 1

#include "merlin_config.h"
// LinearAlgebra
#include "TLAS/LinearAlgebra.h"
// Space3D
#include "EuclideanGeometry/Space3D.h"
// RotationType
#include "EuclideanGeometry/RotationType.h"
// Rot3Drep
#include "EuclideanGeometry/Rot3Drep.h"

//	A general 3 dimensional rotation. Rotation3D objects can
//	be multiplied together, and can also act on points or
//	vectors in space.  The rotations can be thought of as a
//	rotation of a set of right-handed coodinate axes.

class Rotation3D
{
public:

    //	Default constructor. Creates an identity (null) rotation.
    Rotation3D ();

    //	Copy constructor.
    Rotation3D (const Rotation3D& R);

    //	Destructor.
    ~Rotation3D ();

    //	Copy assignment.
    const Rotation3D& operator = (const Rotation3D& right);

    //	Rotate the specified point and return the result.
    Point3D operator () (const Point3D& x) const;

    //	Rotate the specified vector and return the result.
    Vector3D operator () (const Vector3D& v) const;

    //	Concatenate this rotation with the argument, returning
    //	the result.
    Rotation3D operator * (const Rotation3D& R) const;

    //	Multiplication form of Point3D rotation.
    Point3D operator * (const Point3D& p) const;

    //	Multiplication form (dot) of Vector3D rotation.
    Vector3D operator * (const Vector3D& v) const;

    //	Return the inverse of this rotation.
    Rotation3D inv () const;

    //	Return the type of rotation. Can be identity, xrot,
    //	yrot, zrot or general.
    RotationType type () const;

    //	Returns true if this is a null rotation.
    bool isIdentity () const;

    //	Returns true if a pure rotation about the x-axis.
    bool isXrot () const;

    //	Returns true if a pure rotation about the y-axis.
    bool isYrot () const;

    //	Returns true if a pure rotation about the z-axis.
    bool isZrot () const;

    //	Static factory method. Constructs an identity (null)
    //	rotation.
    static Rotation3D identity ();

    //	Static factory method. Constructs a pure x rotation of
    //	angle radians.
    static Rotation3D rotationX (double angle);

    //	Static factory method. Constructs a pure y rotation of
    //	angle radians.
    static Rotation3D rotationY (double angle);

    //	Static factory method. Constructs a pure z rotation of
    //	angle radians.
    static Rotation3D rotationZ (double angle);

    //	Special functions which rotate this transformation by
    //	180 degrees about the specified axis. If R is the 180
    //	degree rotation, and *this is T, then T <- R.T.R. These
    //	functions are supplied for efficiency reasons.
    Rotation3D rotXbyPI () const;
    Rotation3D rotYbyPI () const;
    Rotation3D rotZbyPI () const;

    //	Returns in m the 3x3 rotation matrix.
    RealMatrix& getMatrix (RealMatrix& m) const;

private:

    Rot3Drep* rep;

    //	Special implementation constructor taking a newly
    //	constructed Rot3Drep object.
    Rotation3D (Rot3Drep* nrep);
};

inline Rotation3D::Rotation3D ()
        : rep(Rot3Drep::identity())
{
    rep->refc++;
}

inline Rotation3D::Rotation3D (const Rotation3D& R)
{
    (rep=R.rep)->refc++;
}

inline Rotation3D::Rotation3D (Rot3Drep* nrep)
{
    (rep=nrep)->refc++;
}

inline Rotation3D::~Rotation3D ()
{
    if(--rep->refc==0) delete rep;
}

inline const Rotation3D& Rotation3D::operator = (const Rotation3D& right)
{
    if(rep!=right.rep) {
        if(--rep->refc==0) delete rep;
        (rep=right.rep)->refc++;
    }
    return *this;
}

inline Point3D Rotation3D::operator () (const Point3D& x) const
{
    return rep->rotate(x);
}

inline Vector3D Rotation3D::operator () (const Vector3D& v) const
{
    return rep->rotate(v);
}

inline Rotation3D Rotation3D::operator * (const Rotation3D& R) const
{
    return Rotation3D(rep->dot(*R.rep));
}

inline Point3D Rotation3D::operator * (const Point3D& p) const
{
    return rep->rotate(p);
}

inline Vector3D Rotation3D::operator * (const Vector3D& v) const
{
    return rep->rotate(v);
}

inline Rotation3D Rotation3D::inv () const
{
    return Rotation3D(rep->inv());
}

inline RotationType Rotation3D::type () const
{
    return rep->type();
}

inline bool Rotation3D::isIdentity () const
{
    return rep->isIdentity();
}

inline bool Rotation3D::isXrot () const
{
    return rep->isXrot();
}

inline bool Rotation3D::isYrot () const
{
    return rep->isYrot();
}

inline bool Rotation3D::isZrot () const
{
    return rep->isZrot();
}

inline Rotation3D Rotation3D::identity ()
{
    return Rotation3D();
}

inline Rotation3D Rotation3D::rotationX (double angle)
{
    return Rotation3D(Rot3Drep::rotationX(angle));
}

inline Rotation3D Rotation3D::rotationY (double angle)
{
    return Rotation3D(Rot3Drep::rotationY(angle));
}

inline Rotation3D Rotation3D::rotationZ (double angle)
{
    return Rotation3D(Rot3Drep::rotationZ(angle));
}

inline Rotation3D Rotation3D::rotXbyPI () const
{
    return Rotation3D(rep->rotXbyPI());
}

inline Rotation3D Rotation3D::rotYbyPI () const
{
    return Rotation3D(rep->rotYbyPI());
}

inline Rotation3D Rotation3D::rotZbyPI () const
{
    return Rotation3D(rep->rotZbyPI());
}

inline RealMatrix& Rotation3D::getMatrix (RealMatrix& m) const
{
    return rep->getMatrix(m);
}

#endif
