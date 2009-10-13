/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (2000)
* 
* file Merlin\BasicTransport\PSvectorTransform3D.h
* last modified 04/04/01 14:41:34
*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
*
* Copyright (c) 2000 by The Merlin Collaboration.  
* ALL RIGHTS RESERVED. 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/

#include "BeamDynamics/SMPTracking/SMPTransform3D.h"
#include "NumericalUtils/MatrixPrinter.h"

namespace {
using namespace SMPTracking;

void RotateMoments(const R2Map& R, PSmoments4D& S)
{
    // Here we assume that the rotations are small, representing
    // typical values of alignment errors. For pure Z rotations the
    // result in exact for any angle. For mixed rotations of pure x and
    // y rotations, the results are not exact, and are only approximately
    // valid for small rotations (angle<<1 radian)

    // centroid (first-order moments)
    const double x = R.r11*S.x()+R.r12*S.y();
    const double xp= R.r11*S.xp()+R.r12*S.yp();
    const double y = R.r21*S.x()+R.r22*S.y();
    const double yp= R.r21*S.xp()+R.r22*S.yp();
    S.x()=x;
    S.xp()=xp;
    S.y()=y;
    S.yp()=yp;

    // second-order moments
    const double s00=R.r11*(R.r11*S(0,0) + R.r12*S(2,0)) + R.r12*(R.r11*S(2,0) + R.r12*S(2,2));
    const double s10=R.r11*(R.r11*S(1,0) + R.r12*S(3,0)) + R.r12*(R.r11*S(2,1) + R.r12*S(3,2));
    const double s11=R.r11*(R.r11*S(1,1) + R.r12*S(3,1)) + R.r12*(R.r11*S(3,1) + R.r12*S(3,3));
    const double s20=R.r11*(R.r21*S(0,0) + R.r22*S(2,0)) + R.r12*(R.r21*S(2,0) + R.r22*S(2,2));
    const double s21=R.r11*(R.r21*S(1,0) + R.r22*S(2,1)) + R.r12*(R.r21*S(3,0) + R.r22*S(3,2));
    const double s22=R.r21*(R.r21*S(0,0) + R.r22*S(2,0)) + R.r22*(R.r21*S(2,0) + R.r22*S(2,2));
    const double s30=R.r11*(R.r21*S(1,0) + R.r22*S(3,0)) + R.r12*(R.r21*S(2,1) + R.r22*S(3,2));
    const double s31=R.r11*(R.r21*S(1,1) + R.r22*S(3,1)) + R.r12*(R.r21*S(3,1) + R.r22*S(3,3));
    const double s32=R.r21*(R.r21*S(1,0) + R.r22*S(3,0)) + R.r22*(R.r21*S(2,1) + R.r22*S(3,2));
    const double s33=R.r21*(R.r21*S(1,1) + R.r22*S(3,1)) + R.r22*(R.r21*S(3,1) + R.r22*S(3,3));
    S(0,0)=s00;
    S(1,0)=s10;
    S(1,1)=s11;
    S(2,0)=s20;
    S(2,1)=s21;
    S(2,2)=s22;
    S(3,0)=s30;
    S(3,1)=s31;
    S(3,2)=s32;
    S(3,3)=s33;
}
};

namespace SMPTracking {

SMPTransform3D::SMPTransform3D (const Transform3D& t)
        : bNoRot(t.R().isIdentity()),nullRotation(false)
{
    if(!bNoRot) {

        RealMatrix RR(3,3);
        t.R().getMatrix(RR);
        //			{MatrixForm(RR,cout);char c; cin.get(c);}
        R.r11=RR(0,0);
        R.r12=RR(0,1);
        R.r21=RR(1,0);
        R.r22=RR(1,1);
        theta_x = RR(0,2);
        theta_y = RR(1,2);

        nullRotation = R.IsIdentity();
    }
    delta_x = t.X().x;
    delta_y = t.X().y;
}


// Apply (approximate) transformation
SliceMacroParticle& SMPTransform3D::Apply(SliceMacroParticle& p) const
{
    p.x()-=delta_x;
    p.y()-=delta_y;

    if(!bNoRot) {
        if(!nullRotation)
            RotateMoments(R,p);
        p.xp()+=theta_x;
        p.yp()+=theta_y;
    }

    return p;
}

}; // end namespace SMPTracking

