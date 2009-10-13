/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:29:36 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
// LatticeFrame
#include "AcceleratorModel/Frames/LatticeFrame.h"

#define VALID_SFRAME(sframe) \
assert(sframe==this || sframe==GLOBAL_FRAME || superFrame!=GLOBAL_FRAME)

Transform3D LatticeFrame::GetFrameTransform (const LatticeFrame* sframe) const
{

    if(sframe==GLOBAL_FRAME)
        sframe = GetGlobalFrame();

    // The frame transformation to itself is the identity
    if(sframe==this)
        return Transform3D();

    if(sframe==superFrame)
        return GetLocalFrameTransform();


    Transform3D t1 = GetPhysicalTransform(sframe);
    double s = GetPosition(sframe);


    Transform3D t2 = sframe->GetGeometryTransform(s);

    return t1*t2.inv();
}

Transform3D LatticeFrame::GetPhysicalTransform (const LatticeFrame* sframe) const
{
    VALID_SFRAME(sframe);

    // The frame transformation to itself is the identity.
    if(sframe==this)
        return Transform3D();

    Transform3D t0 = GetLocalFrameTransform();

    if(superFrame!=GLOBAL_FRAME) {
        Transform3D t1 = superFrame->GetGeometryTransform(0,s_0);
        Transform3D t2 = superFrame->GetPhysicalTransform(sframe);
        t0=t0*t1*t2;
    }

    return t0;
}

double LatticeFrame::GetPosition (const LatticeFrame* sframe) const
{
    VALID_SFRAME(sframe);
    return sframe==superFrame ? s_0 : s_0 + superFrame->GetPosition(sframe);
}

AcceleratorGeometry::Extent LatticeFrame::GetGeometryExtent (const LatticeFrame* sframe) const
{
    AcceleratorGeometry::Extent extent = GetLocalGeometryExtent();
    double s = GetPosition(sframe);
    extent.first += s;
    extent.second += s;
    return extent;
}

Transform3D LatticeFrame::GetBoundaryPlaneTransform (BoundaryPlane p) const
{
    // we ignore any global transformations of the top-level frame
    if(superFrame==GLOBAL_FRAME)
        return Transform3D();

    Transform3D t0=LocalBoundaryPlaneTransform(p);
    if(superFrame!=0 && superFrame->IsBoundaryPlane(p,this))
        t0=t0*(superFrame->GetBoundaryPlaneTransform(p));
    return t0;
}
/****
void LatticeFrame::Translate (double dx, double dy, double dz)
{
    _TRNSFM(translation(dx,dy,dz));
}

void LatticeFrame::RotateX (double angle)
{
    _TRNSFM(rotationX(angle));
}

void LatticeFrame::RotateY (double angle)
{
    _TRNSFM(rotationY(angle));
}

void LatticeFrame::RotateZ (double angle)
{
    _TRNSFM(rotationZ(angle));
}
***/

void LatticeFrame::ApplyLocalFrameTransform (const Transform3D& t)
{
    if(!local_T)
        local_T = new Transform3D(t);
    else
        (*local_T)*=t;

    Invalidate();
}

void LatticeFrame::SetLocalFrameTransform (const Transform3D& t)
{
    if(!local_T)
        local_T = new Transform3D(t);
    else
        (*local_T)=t;

    Invalidate();
}

void LatticeFrame::ClearLocalFrameTransform ()
{
    ClearTransform();
}

LatticeFrame* LatticeFrame::GetGlobalFrame () const
{
    return superFrame==0 ? const_cast<LatticeFrame*>(this) : superFrame->GetGlobalFrame();
}

Transform3D LatticeFrame::LocalBoundaryPlaneTransform (BoundaryPlane p) const
{
    Transform3D t0 = GetLocalFrameTransform();
    if(!t0.isIdentity()) {
        Transform3D t = GetGeometryTransform(p);
        return t*t0*t.inv();
    }
    return t0;
}

void LatticeFrame::Traverse(FrameTraverser &ft)
{
    ft.ActOn(this);
}
