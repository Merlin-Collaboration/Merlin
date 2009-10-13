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

// memcpy
#include <cstring>
// blas template routines
#include "TLAS/LinearAlgebra.h"
// rotation matrix construction
#include "EuclideanGeometry/RotationMatrix.h"
// GeneralRotation
#include "EuclideanGeometry/GeneralRotation.h"

GeneralRotation::GeneralRotation (const GeneralRotation& gr)
        :m(gr.m)
{}

GeneralRotation::GeneralRotation (const RotationX& rx, const RotationY& ry)
        :m(3,3)
{
    m(0,0) = ry.cosine();
    m(0,1) = 0;
    m(0,2) = -ry.sine();

    m(1,0) = rx.sine()*ry.sine();
    m(1,1) = rx.cosine();
    m(1,2) = rx.sine()*ry.cosine();

    m(2,0) = rx.cosine()*ry.sine();
    m(2,1) = -rx.sine();
    m(2,2) = rx.cosine()*ry.cosine();
}

GeneralRotation::GeneralRotation (const RotationX& rx, const RotationZ& rz)
        :m(3,3)
{
    m(0,0) = rz.cosine();
    m(0,1) = rz.sine();
    m(0,2) = 0;

    m(1,0) = -rx.cosine()*rz.sine();
    m(1,1) = rx.cosine()*rz.cosine();
    m(1,2) = rx.sine();

    m(2,0) = rx.sine()*rz.sine();
    m(2,1) = -rz.cosine()*rx.sine();
    m(2,2) = rx.cosine();
}

GeneralRotation::GeneralRotation (const RotationY& ry, const RotationX& rx)
        :m(3,3)
{
    m(0,0) = ry.cosine();
    m(0,1) = rx.sine()*ry.sine();
    m(0,2) = -rx.cosine()*ry.sine();

    m(1,0) = 0;
    m(1,1) = rx.cosine();
    m(1,2) = rx.sine();

    m(2,0) = ry.sine();
    m(2,1) = -ry.cosine()*rx.sine();
    m(2,2) = rx.cosine()*ry.cosine();
}

GeneralRotation::GeneralRotation (const RotationY& ry, const RotationZ& rz)
        :m(3,3)
{
    m(0,0) = ry.cosine()*rz.cosine();
    m(0,1) = ry.cosine()*rz.sine();
    m(0,2) = -ry.sine();

    m(1,0) = -rz.sine();
    m(1,1) = rz.cosine();
    m(1,2) = 0;

    m(2,0) = rz.cosine()*ry.sine();
    m(2,1) = ry.sine()*rz.sine();
    m(2,2) = ry.cosine();
}

GeneralRotation::GeneralRotation (const RotationZ& rz, const RotationX& rx)
        :m(3,3)
{
    m(0,0) = rz.cosine();
    m(0,1) = rx.cosine()*rz.sine();
    m(0,2) = rx.sine()*rz.sine();

    m(1,0) = -rz.sine();
    m(1,1) = rx.cosine()*rz.cosine();
    m(1,2) = rz.cosine()*rx.sine();

    m(2,0) = 0;
    m(2,1) = -rx.sine();
    m(2,2) = rx.cosine();
}

GeneralRotation::GeneralRotation (const RotationZ& rz, const RotationY& ry)
        :m(3,3)
{
    m(0,0) = ry.cosine()*rz.cosine();
    m(0,1) = rz.sine();
    m(0,2) = -rz.cosine()*ry.sine();

    m(1,0) = -ry.cosine()*rz.sine();
    m(1,1) = rz.cosine();
    m(1,2) = ry.sine()*rz.sine();

    m(2,0) = ry.sine();
    m(2,1) = 0;
    m(2,2) = ry.cosine();
}

GeneralRotation::GeneralRotation (const GeneralRotation& r, const RotationX& rx)
        : m(r.m*DECL_XROT(rx.cosine(),rx.sine()))
{}

GeneralRotation::GeneralRotation (const GeneralRotation& r, const RotationY& ry)
        : m(r.m*DECL_YROT(ry.cosine(),ry.sine()))
{}

GeneralRotation::GeneralRotation (const GeneralRotation& r, const RotationZ& rz)
        : m(r.m*DECL_ZROT(rz.cosine(),rz.sine()))
{}

Rot3Drep* GeneralRotation::inv () const
{
    GeneralRotation *newrot = new GeneralRotation;
    for(int i=0; i<3; i++)
        for(int j=0; j<=i; j++) {
            if(i==j)
                newrot->m(i,i) = m(i,i);
            else {
                newrot->m(i,j) = m(j,i);
                newrot->m(j,i) = m(i,j);
            }
        }
    return newrot;
}

Point3D GeneralRotation::rotate (const Point3D& p) const
{
    double x = m(0,0)*p.x+m(0,1)*p.y+m(0,2)*p.z;
    double y = m(1,0)*p.x+m(1,1)*p.y+m(1,2)*p.z;
    double z = m(2,0)*p.x+m(2,1)*p.y+m(2,2)*p.z;
    return Point3D(x,y,z);
}

Vector3D GeneralRotation::rotate (const Vector3D& v) const
{
    double x = m(0,0)*v.x+m(0,1)*v.y+m(0,2)*v.z;
    double y = m(1,0)*v.x+m(1,1)*v.y+m(1,2)*v.z;
    double z = m(2,0)*v.x+m(2,1)*v.y+m(2,2)*v.z;
    return Vector3D(x,y,z);
}

Rot3Drep* GeneralRotation::dot (const Rot3Drep& r) const
{
    return r.dotBy(*this);
}

Rot3Drep* GeneralRotation::dotBy (const RotationX& rx) const
{
    return new GeneralRotation(DECL_XROT(rx.cosine(),rx.sine())*m);
}

Rot3Drep* GeneralRotation::dotBy (const RotationY& ry) const
{
    return new GeneralRotation(DECL_YROT(ry.cosine(),ry.sine())*m);
}

Rot3Drep* GeneralRotation::dotBy (const RotationZ& rz) const
{
    return new GeneralRotation(DECL_ZROT(rz.cosine(),rz.sine())*m);
}

Rot3Drep* GeneralRotation::dotBy (const GeneralRotation& r) const
{
    return new GeneralRotation((r.m)*m);
}

RotationType GeneralRotation::type () const
{
    return mixed;
}

bool GeneralRotation::isIdentity () const
{
    return m(0,0)==1.0 && m(1,1)==1.0 && m(2,2)==1.0;
}

bool GeneralRotation::isXrot () const
{
    return m(0,0)==1.0 && m(0,1)==0.0 && m(0,2)==0.0;
}

bool GeneralRotation::isYrot () const
{
    return m(1,0)==0.0 && m(1,1)==1.0 && m(1,2)==0.0;
}

bool GeneralRotation::isZrot () const
{
    return m(2,0)==0.0 && m(2,1)==0.0 && m(2,2)==1.0;
}

Rot3Drep* GeneralRotation::rotXbyPI () const
{
    RealMatrix R(m);
    R(0,1)*=-1;
    R(0,2)*=-1;
    R(1,0)*=-1;
    R(2,0)*=-1;
    return new GeneralRotation(R);
}

Rot3Drep* GeneralRotation::rotYbyPI () const
{
    RealMatrix R(m);
    R(0,1)*=-1;
    R(1,0)*=-1;
    R(1,2)*=-1;
    R(2,1)*=-1;
    return new GeneralRotation(R);
}

Rot3Drep* GeneralRotation::rotZbyPI () const
{
    RealMatrix R(m);
    R(0,2)*=-1;
    R(1,2)*=-1;
    R(2,0)*=-1;
    R(2,1)*=-1;
    return new GeneralRotation(R);
}

RealMatrix& GeneralRotation::getMatrix (RealMatrix& m) const
{
    return m.copy(this->m);
}

