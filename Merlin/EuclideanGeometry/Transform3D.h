


#ifndef Transform3D_h
#define Transform3D_h 1

#include "merlin_config.h"


// Transform2D
#include "EuclideanGeometry/Transform2D.h"
// Space3D
#include "EuclideanGeometry/Space3D.h"
// Rotation3D
#include "EuclideanGeometry/Rotation3D.h"



//	Represents a 3-dimensional Euclidean transformation. A
//	Euclidean transformation is characterised by a
//	translation followed by a rotation. The translation
//	vector and rotation can be thought of as the origin and
//	orientation respectively of a set of mutually
//	perpendicular axes. EuclideanTransform objects can act
//	on Point3D and Vector3D objects, effectively
//	transforming them into a new reference frame.
//	EuclideanTransforms can also be composed using the
//	multiplication operator *, i.e. T3 = T1*T2 (note that
//	multiplication is not commutitive).


class Transform3D
{
public:
    //	Default constructor. Creates an indentity (null)
    //	transformation.
    Transform3D ();

    //	Constructor taking explicitely a translation and a
    //	rotation.
    Transform3D (const Point3D& X, const Rotation3D& R);

    //	Construct a general 3D transformation from a Transform2D
    //	object.
    explicit Transform3D (const Transform2D& t2d);


    //	Returns the translation associated with this
    //	transformation.
    const Point3D& X () const;

    //	Returns the rotation associated with this transformation.
    const Rotation3D& R () const;

    //	Transform the supplied Point3D and return the result.
    //	Points are translated and rotated.
    Point3D operator () (const Point3D& p) const;

    //	Transform the supplied Vector3D and return the result.
    //	Vectors are rotated only.
    Vector3D operator () (const Vector3D& v) const;

    //	Multiplication operator form of transformation
    //	operationon a Point3D.
    Point3D operator * (const Point3D& p) const;

    //	Multiplication operator form of transformation operation
    //	on a Vector3D.
    Vector3D operator * (const Vector3D& v) const;

    //	Composition (non-commutative multiplication) operator.
    Transform3D operator * (const Transform3D& t) const;

    //	Composes this transformation with t. If this
    //	transformation is t0, then t0*=t is mathematically
    //	equivalent to t0<-t*t0.
    const Transform3D& operator *= (const Transform3D& t);

    //	Return the inverse of this transformation.
    Transform3D inv () const;

    //	Invert this transformation.
    const Transform3D& invert ();

    //	Special functions which rotate this transformation by
    //	180 degrees about the specified axis. If R is the 180
    //	degree rotation, and *this is T, then T <- R.T.R. These
    //	functions are supplied for efficiency reasons.
    void rotXbyPI ();

    void rotYbyPI ();

    void rotZbyPI ();

    //	Returns true if this is the identity transformation.
    bool isIdentity () const;

    //	Static (factory) method for construcing a pure
    //	x-rotation transformation.
    static Transform3D rotationX (double angle);

    //	Static (factory) method for construcing a pure
    //	y-rotation transformation.
    static Transform3D rotationY (double angle);

    //	Static (factory) method for construcing a pure
    //	z-rotation transformation.
    static Transform3D rotationZ (double angle);

    //	Static (factory) method for construcing a pure
    //	translation transformation.
    static Transform3D translation (double dx, double dy, double dz);

protected:
private:
    // Data Members for Associations

    //	The translation associated with this transformation.
    Point3D x0;

    //	The rotation associated with this transformation.
    Rotation3D r;

private:
};

// Class Transform3D

inline Transform3D::Transform3D ()
        :x0(0,0,0),r(Rotation3D::identity())
{
}

inline Transform3D::Transform3D (const Point3D& X, const Rotation3D& R)
        :x0(X),r(R)
{
}


inline Transform3D::Transform3D (const Transform2D& t2d)
        :x0(t2d.translation().x,t2d.translation().y,0),r(Rotation3D::rotationZ(t2d.rotationAngle()))
{
}



inline const Point3D& Transform3D::X () const
{
    return x0;
}

inline const Rotation3D& Transform3D::R () const
{
    return r;
}

inline Point3D Transform3D::operator () (const Point3D& p) const
{
    return r(p-x0);
}

inline Vector3D Transform3D::operator () (const Vector3D& v) const
{
    return r(v);
}

inline Point3D Transform3D::operator * (const Point3D& p) const
{
    return r(p-x0);
}

inline Vector3D Transform3D::operator * (const Vector3D& v) const
{
    return r(v);
}

inline Transform3D Transform3D::operator * (const Transform3D& t) const
{
    Transform3D ret(t);
    return ret*=(*this);
}

inline void Transform3D::rotXbyPI ()
{
    x0.y*=-1;
    x0.z*=-1;
    r = r.rotXbyPI();
}

inline void Transform3D::rotYbyPI ()
{
    x0.x*=-1;
    x0.z*=-1;
    r = r.rotYbyPI();
}

inline void Transform3D::rotZbyPI ()
{
    x0.x*=-1;
    x0.y*=-1;
    r=r.rotZbyPI();
}

inline bool Transform3D::isIdentity () const
{
    return x0.isZero() && r.isIdentity();
}



#endif
