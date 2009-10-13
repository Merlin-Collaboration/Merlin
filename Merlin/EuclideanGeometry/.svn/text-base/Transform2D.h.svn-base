


#ifndef Transform2D_h
#define Transform2D_h 1

#include "merlin_config.h"

#include <cmath>
#include <iostream>

// Space2D
#include "EuclideanGeometry/Space2D.h"

using std::ostream;
using std::istream;


//	Represents a simpy X-Y rotoation.


class Rotation2D
{
public:
    //	Constructor taking the rotation angle (in radians).
    explicit Rotation2D (double angle = 0)
            : cosa(cos(angle)),sina(sin(angle))
    {
    }


    //	Returns the rotation angle in radian.
    double rotation () const
    {
        return atan2(sina,cosa);
    }

    //	Set the rotation angle.
    void setRotation (double phi)
    {
        cosa=cos(phi);
        sina=sin(phi);
    }

    //	Returns sin(rotation).
    double sinphi () const
    {
        return sina;
    }

    //	Returns cos(rotation).
    double cosphi () const
    {
        return cosa;
    }

    //	Returns true if the rotation angle is zero.
    bool isIdentity () const
    {
        return sina == 0;
    }

    //	Premultipliy *this by R.
    const Rotation2D& operator *= (const Rotation2D& R)
    {
        const double c = R.cosa*cosa-R.sina*sina;
        const double s = R.sina*cosa+R.cosa*sina;
        cosa=c;
        sina=s;
        return *this;
    }

    //	Rotation multiplication.
    friend Rotation2D operator * (const Rotation2D& R2, const Rotation2D& R1)
    {
        return Rotation2D(R2.cosa*R1.cosa-R2.sina*R1.sina,
                          R2.sina*R1.cosa+R2.cosa*R1.sina);
    }

    //	Return the inverse of *this.
    Rotation2D inv () const
    {
        return Rotation2D(cosa,-sina);
    }

    //	Invert (in-place) *this.
    const Rotation2D& invert ()
    {
        sina=-sina;
        return *this;
    }

    //	Transformation operations.
    Point2D operator () (const Point2D& p) const
    {
        return Point2D(cosa*p.x-sina*p.y,sina*p.x+cosa*p.y);
    }

    Vector2D operator () (const Vector2D& v) const
    {
        return Vector2D(cosa*v.x-sina*v.y,sina*v.x+cosa*v.y);
    }

    friend ostream& operator << (ostream& os, const Rotation2D& rot)
    {
        os<<rot.cosa;
        os<<rot.sina;
        return os;
    }

    friend istream& operator >> (istream& is, Rotation2D& rot)
    {
        is>>rot.cosa;
        is>>rot.sina;
        return is;
    }

protected:
private:
    //	Private constructor.
    Rotation2D (double c, double s)
            : cosa(c),sina(s)
    {
    }

    // Data Members for Class Attributes

    //	cosine of the rotation angle.
    double cosa;

    //	sin of the rotation angle.
    double sina;

private:
};

//	A special sub-set of a full 3D Euclidean transformation,
//	which represents a rotation in the Z-plane followed by
//	an X and Y displacement.


class Transform2D
{
public:
    //	Constructor.
    Transform2D (double dx = 0, double dy = 0, double dphi = 0)
            : r(dphi), x0(dx,dy)
    {
    }

    explicit Transform2D (const Rotation2D& rr, const Point2D& p = Point2D(0,0))
            : r(rr),x0(p)
    {
    }

    explicit Transform2D (const Point2D& p)
            : r(0),x0(p)
    {
    }


    //	Multiplication (concatanation) operator.
    friend Transform2D operator * (const Transform2D& t2, const Transform2D& t1)
    {
        // t3=t2*t1;
        return Transform2D(t2.r*t1.r,t2.r(t1.x0)+t2.x0);
    }

    //	t1*=t2 facilitates pre-multiplication of t1 by t2, i.e.
    //	t1->t2*t1.
    const Transform2D& operator *= (const Transform2D& t)
    {
        return *this = t*(*this);
    }

    //	Invert the transformation and return the result.
    Transform2D inv () const
    {
        Transform2D ti(*this);
        return ti.invert();
    }

    //	Invert (in-place) this transformation.
    const Transform2D& invert ()
    {
        r.invert();
        x0 = -r(x0);
        return *this;
    }

    //	Returns the (x,y) translation as a Point2D.
    const Point2D& translation () const
    {
        return x0;
    }

    //	Returns the rotation object.
    const Rotation2D& rotation () const
    {
        return r;
    }

    //	Returns the rotation angle in radian.
    double rotationAngle () const
    {
        return r.rotation();
    }

    //	Returns sin(rotation).
    double sinphi () const
    {
        return r.sinphi();
    }

    //	Returns cos(rotation).
    double cosphi () const
    {
        return r.cosphi();
    }

    //	Returns true if this transform is the identity.
    bool isIdentity () const
    {
        return x0.isZero() && r.isIdentity();
    }

    //	Transformation operations.
    Point2D operator () (const Point2D& p) const
    {
        return r(p)-x0;
    }

    Vector2D operator () (const Vector2D& v) const
    {
        return r(v);
    }

    //	Modifiers.
    void setTranslationX (double x)
    {
        x0.x=x;
    }

    void setTranslationY (double y)
    {
        x0.y=y;
    }

    void setRotation (double phi)
    {
        r.setRotation(phi);
    }

    friend ostream& operator << (ostream& os, const Transform2D& t)
    {
        os<<t.x0.x;
        os<<t.x0.y;
        os<<t.r;
        return os;
    }

    friend istream& operator >> (istream& is, Transform2D& t)
    {
        is>>t.x0.x;
        is>>t.x0.y;
        is>>t.r;
        return is;
    }

protected:
private:
    // Data Members for Associations

    Rotation2D r;

    Point2D x0;

private:
};

// Class Rotation2D




// Class Transform2D






#endif
