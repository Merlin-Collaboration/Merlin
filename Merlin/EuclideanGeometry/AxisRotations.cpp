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

#include "EuclideanGeometry/IdentityRotation.h"
#include "EuclideanGeometry/RotationMatrix.h"
// GeneralRotation
#include "EuclideanGeometry/GeneralRotation.h"
// AxisRotations
#include "EuclideanGeometry/AxisRotations.h"

double PureAxisRotation::angle () const
{
    return atan2(sinphi,cosphi);
}

Rot3Drep* RotationX::inv () const
{
    return new RotationX(cosine(),-sine());
}

Point3D RotationX::rotate (const Point3D& p) const
{
    return Point3D(p.x,p.y*cosine()+p.z*sine(),-p.y*sine()+p.z*cosine());
}

Vector3D RotationX::rotate (const Vector3D& v) const
{
    return Vector3D(v.x,v.y*cosine()+v.z*sine(),-v.y*sine()+v.z*cosine());
}

Rot3Drep* RotationX::dot (const Rot3Drep& r) const
{
    return r.dotBy(*this);
}

Rot3Drep* RotationX::dotBy (const RotationX& rx) const
{
    double ncos = cosine()*rx.cosine()-sine()*rx.sine();
    if(ncos==1)
        return new IdentityRotation();
    else {
        double nsin = sine()*rx.cosine()+cosine()*rx.sine();
        return new RotationX(ncos,nsin);
    }
}

Rot3Drep* RotationX::dotBy (const RotationY& ry) const
{
    return new GeneralRotation(ry,*this);
}

Rot3Drep* RotationX::dotBy (const RotationZ& rz) const
{
    return new GeneralRotation(rz,*this);
}

Rot3Drep* RotationX::dotBy (const GeneralRotation& r) const
{
    return new GeneralRotation(r,*this);
}

RotationType RotationX::type () const
{
    return xrot;
}

bool RotationX::isXrot () const
{
    return true;
}

Rot3Drep* RotationX::rotXbyPI () const
{
    return const_cast<RotationX*>(this);
}

Rot3Drep* RotationX::rotYbyPI () const
{
    return new RotationX(cosine(),-sine());
}

Rot3Drep* RotationX::rotZbyPI () const
{
    return new RotationX(cosine(),-sine());
}

RealMatrix& RotationX::getMatrix (RealMatrix& m) const
{
    return m.copy(DECL_XROT(cosine(),sine()));
}

// Class RotationY


Rot3Drep* RotationY::inv () const
{
    return new RotationY(cosine(),-sine());
}

Point3D RotationY::rotate (const Point3D& p) const
{
    return Point3D(p.x*cosine()-p.z*sine(),p.y,p.z*cosine()+p.x*sine());
}

Vector3D RotationY::rotate (const Vector3D& v) const
{
    return Vector3D(v.x*cosine()-v.z*sine(),v.y,v.z*cosine()+v.x*sine());
}

Rot3Drep* RotationY::dot (const Rot3Drep& r) const
{
    return r.dotBy(*this);
}

Rot3Drep* RotationY::dotBy (const RotationX& rx) const
{
    return new GeneralRotation(rx,*this);
}

Rot3Drep* RotationY::dotBy (const RotationY& ry) const
{
    double ncos = cosine()*ry.cosine()-sine()*ry.sine();
    if(ncos==1)
        return new IdentityRotation();
    else {
        double nsin = sine()*ry.cosine()+cosine()*ry.sine();
        return new RotationY(ncos,nsin);
    }
}

Rot3Drep* RotationY::dotBy (const RotationZ& rz) const
{
    return new GeneralRotation(rz,*this);
}

Rot3Drep* RotationY::dotBy (const GeneralRotation& r) const
{
    return new GeneralRotation(r,*this);
}

RotationType RotationY::type () const
{
    return yrot;
}

bool RotationY::isYrot () const
{
    return true;
}

Rot3Drep* RotationY::rotXbyPI () const
{
    return new RotationY(cosine(),-sine());
}

Rot3Drep* RotationY::rotYbyPI () const
{
    return const_cast<RotationY*>(this);
}

Rot3Drep* RotationY::rotZbyPI () const
{
    return new RotationY(cosine(),-sine());
}

RealMatrix& RotationY::getMatrix (RealMatrix& m) const
{
    return m.copy(DECL_YROT(cosine(),sine()));
}

// Class RotationZ


Rot3Drep* RotationZ::inv () const
{
    return new RotationZ(cosine(),-sine());
}

Point3D RotationZ::rotate (const Point3D& p) const
{
    return Point3D(p.x*cosine()+p.y*sine(),p.y*cosine()-p.x*sine(),p.z);
}

Vector3D RotationZ::rotate (const Vector3D& v) const
{
    return Vector3D(v.x*cosine()+v.y*sine(),v.y*cosine()-v.x*sine(),v.z);
}

Rot3Drep* RotationZ::dot (const Rot3Drep& r) const
{
    return r.dotBy(*this);
}

Rot3Drep* RotationZ::dotBy (const RotationX& rx) const
{
    return new GeneralRotation(rx,*this);
}

Rot3Drep* RotationZ::dotBy (const RotationY& ry) const
{
    return new GeneralRotation(ry,*this);
}

Rot3Drep* RotationZ::dotBy (const RotationZ& rz) const
{
    double ncos = cosine()*rz.cosine()-sine()*rz.sine();
    if(ncos==1)
        return new IdentityRotation();
    else {
        double nsin = sine()*rz.cosine()+cosine()*rz.sine();
        return new RotationZ(ncos,nsin);
    }
}

Rot3Drep* RotationZ::dotBy (const GeneralRotation& r) const
{
    return new GeneralRotation(r,*this);
}

RotationType RotationZ::type () const
{
    return zrot;
}

bool RotationZ::isZrot () const
{
    return true;
}

Rot3Drep* RotationZ::rotXbyPI () const
{
    return new RotationZ(cosine(),-sine());
}

Rot3Drep* RotationZ::rotYbyPI () const
{
    return new RotationZ(cosine(),-sine());
}

Rot3Drep* RotationZ::rotZbyPI () const
{
    return const_cast<RotationZ*>(this);
}

RealMatrix& RotationZ::getMatrix (RealMatrix& m) const
{
    return m.copy(DECL_ZROT(cosine(),sine()));
}

