/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
// PSTypes
#include "BeamModel/PSTypes.h"
// MatrixMaps
#include "BasicTransport/MatrixMaps.h"

namespace {

PSmoments& TransformMoments(PSmoments& S, const RealMatrix& R)
{
    PSmoments result; // initialised to zero
    register int i,j,k,l;
    const int n=R.ncols();

    for(i=0; i<n; i++) {
        for(j=0; j<=i; j++) {
            double sij=0.0;
            for(k=0; k<n; k++)
                for(l=0; l<n; l++)
                    sij+=R(i,k)*R(j,l)*S(k,l);

            result(i,j)=sij;
        }
        for(j=0;j<n;j++)
            result[i]+=R(i,j)*S[j]; // centroid tracking
    }

    //		for(i=0;i<n;i++)
    //			for(j=0;j<=i;j++)
    //				S(i,j)=result(i,j);

    return S=result;
}

void RescaleMoments(double P0, double P1, PSmoments& S)
{
    // rescale the centroid
    double ratio = P0/P1;
    S[5]=ratio*(1+S[5])-1.0;

    // rescale momentum correlations
    for(int i=0; i<4; i++)
        S(i,5)*=ratio;

    // rescale energy spread
    S(5,5)*=(ratio*ratio);
}

}; // end of anonymous namespace


PSvector& RMtrx::Apply (PSvector& x) const
{
    PSvector x1(0);
    tblas2::tgemv(false,1.0,R,x,0.0,x1);
    for(Subscript i=0;i<R.nrows();i++) // only copy back that which changed
        x[i]=x1[i];
    return x;
}

PSvector& RMtrx::Apply (PSvector& x, double p0) const
{
    if(EnergyIndependent())
        return Apply(x);

    double dp=x.dp();
    x.dp() = scaledp(p0,dp);
    Apply(x);
    x.dp()=dp;
    return x;
}

PSvectorArray& RMtrx::Apply (PSvectorArray& xa) const
{
    for(PSvectorArray::iterator p=xa.begin(); p!=xa.end(); p++)
        Apply(*p);
    return xa;
}

PSvectorArray& RMtrx::Apply (PSvectorArray& xa, double p0) const
{
    if(EnergyIndependent())
        return Apply(xa);

    for(PSvectorArray::iterator p=xa.begin(); p!=xa.end(); p++)
        Apply(*p,p0);
    return xa;
}

PSmoments& RMtrx::Apply (PSmoments& sigma) const
{
    return TransformMoments(sigma,R);
}

PSmoments& RMtrx::Apply (PSmoments& sigma, double p0) const
{
    if(EnergyIndependent())
        return Apply(sigma);

    RescaleMoments(p0,Pref,sigma);
    TransformMoments(sigma,R);
    RescaleMoments(Pref,p0,sigma);
    return sigma;
}

PSmomentsArray& RMtrx::Apply (PSmomentsArray& sigmaArray) const
{
    for(PSmomentsArray::iterator p=sigmaArray.begin(); p!=sigmaArray.end(); p++)
        Apply(*p);
    return sigmaArray;
}

PSmomentsArray& RMtrx::Apply (PSmomentsArray& sigmaArray, double p0) const
{
    if(EnergyIndependent())
        return Apply(sigmaArray);

    for(PSmomentsArray::iterator p=sigmaArray.begin(); p!=sigmaArray.end(); p++)
        Apply(*p,p0);
    return sigmaArray;
}

RMtrx& RMtrx::Invert ()
{
    // TODO:
    assert(false);
    return *this;
}

RMtrx& RMtrx::operator *= (const RMtrx& rhs)
{
    R = (rhs.R)*R;
    return *this;
}

PSvector& RdpMtrx::Apply (PSvector& x) const
{
    PSvector x1(0.0);
    int n=R.nrows();
    register int i,j;

    for(i=0;i<n;i++) {
        x1[i]=0;
        for(j=0; j<n; j++)
            x1[i]+=(R(i,j)+T(i,j)*x.dp())*x[j];
    }
    for(i=0;i<n;i++)
        x[i]=x1[i];

    return x;
}

PSvector& RdpMtrx::Apply (PSvector& x, double p0) const
{
    if(EnergyIndependent())
        return Apply(x);

    const double dp=x.dp();
    x.dp()=scaledp(p0,dp);
    Apply(x);
    x.dp()=dp;
    return x;
}

PSvectorArray& RdpMtrx::Apply (PSvectorArray& xa) const
{
    for(PSvectorArray::iterator p=xa.begin(); p!=xa.end(); p++)
        Apply(*p);
    return xa;
}

PSvectorArray& RdpMtrx::Apply (PSvectorArray& xa, double p0) const
{
    if(EnergyIndependent())
        return Apply(xa);

    for(PSvectorArray::iterator p=xa.begin(); p!=xa.end(); p++)
        Apply(*p,p0);
    return xa;
}

PSmoments& RdpMtrx::Apply (PSmoments& sigma) const
{
    RealMatrix R1=R;
    R1+=T*sigma[5];
    return TransformMoments(sigma,R1);
}

PSmoments& RdpMtrx::Apply (PSmoments& sigma, double p0) const
{
    if(EnergyIndependent())
        return Apply(sigma);

    RescaleMoments(p0,Pref,sigma);
    Apply(sigma);
    RescaleMoments(Pref,p0,sigma);
    return sigma;
}

PSmomentsArray& RdpMtrx::Apply (PSmomentsArray& sigmaArray) const
{
    for(PSmomentsArray::iterator p=sigmaArray.begin(); p!=sigmaArray.end(); p++)
        Apply(*p);
    return sigmaArray;
}

PSmomentsArray& RdpMtrx::Apply (PSmomentsArray& sigmaArray, double p0) const
{
    if(EnergyIndependent())
        return Apply(sigmaArray);

    for(PSmomentsArray::iterator p=sigmaArray.begin(); p!=sigmaArray.end(); p++)
        Apply(*p,p0);
    return sigmaArray;
}

LinMtrxBase::LinMtrxBase (const RealMatrix& RR, double P0)
        : R(RR),Pref(P0)
{
#ifndef NDEBUG
    size_t n=R.nrows()/2;
    assert((R.ncols()==2*n)&&(n==1||n==2||n==3));
#endif
}


