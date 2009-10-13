



// Transform3D
#include "EuclideanGeometry/Transform3D.h"



// Class Transform3D





const Transform3D& Transform3D::operator *= (const Transform3D& t)
{
    // this operation represents t*this
    x0+=r.inv()*t.x0;
    r = t.r*r;
    return *this;
}

Transform3D Transform3D::inv () const
{
    return Transform3D(-(r*x0),r.inv());
}

const Transform3D& Transform3D::invert ()
{
    x0 = -(r*x0);
    r = r.inv();
    return *this;
}

Transform3D Transform3D::rotationX (double angle)
{
    return Transform3D(Point3D(0,0,0),Rotation3D::rotationX(angle));
}

Transform3D Transform3D::rotationY (double angle)
{
    return Transform3D(Point3D(0,0,0),Rotation3D::rotationY(angle));
}

Transform3D Transform3D::rotationZ (double angle)
{
    return Transform3D(Point3D(0,0,0),Rotation3D::rotationZ(angle));
}

Transform3D Transform3D::translation (double dx, double dy, double dz)
{
    return Transform3D(Point3D(dx,dy,dz),Rotation3D::identity());
}

