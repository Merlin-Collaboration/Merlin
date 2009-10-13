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
// $Revision: 1.6 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/MatrixPrinter.h"
#include "NumericalUtils/utils.h"
#include "BasicTransport/TransportMatrix.h"

namespace {

// utility functions for matrix calculations

inline RealMatrix& InitR(RealMatrix& R)
{
    R=IdentityMatrix(R.nrows());
    return R;
}

// Note that the following routines assume that
// R is the identity on entrance.

void SROT(double cosphi, double sinphi, RealMatrix& R)
{
    //assert(fequal(cosphi*cosphi,sinphi*sinphi) && R.nrows()>3);
    R(0,0)=R(1,1)=R(2,2)=R(3,3)=cosphi;
    R(0,2)=R(1,3)= -sinphi;
    R(2,0)=R(3,1)= sinphi;
}

void SBR(double length, double h, double Kx, RealMatrix& R)
{
    // check for zero length element
    if(!length)
        return;

    assert(R.nrows()>3);

    double Ky = -Kx;
    Kx += h*h;

    // horizontal plane

    if(Kx==0)
        // simply horizontal drift space
        R(0,1) = length;
    else if(Kx>0) {
        // horizontally focusing magnet
        double k = sqrt(Kx);
        double phi = k*length;
        double sinPhi = sin(phi);

        R(0,0) = R(1,1) = cos(phi);
        R(0,1) = sinPhi/k;
        R(1,0) = -k*sinPhi;
    }
    else {
        // horizontally defocusing magnet
        double k = sqrt(-Kx);
        double phi = k*length;
        double sinhPhi = sinh(phi);

        R(0,0) = R(1,1) = cosh(phi);
        R(0,1) = sinhPhi/k;
        R(1,0) = k*sinhPhi;
    }

    // vectircal plane

    if(Ky==0)
        // simply horizontal drift space
        R(2,3) = length;
    else if(Ky>0) {
        // vertically focusing magnet
        double k = sqrt(Ky);
        double phi = k*length;
        double sinPhi = sin(phi);

        R(2,2) = R(3,3) = cos(phi);
        R(2,3) = sinPhi/k;
        R(3,2) = -k*sinPhi;
    }
    else {
        // vertically defocusing magnet
        double k = sqrt(-Ky);
        double phi = k*length;
        double sinhPhi = sinh(phi);

        R(2,2) = R(3,3) = cosh(phi);
        R(2,3) = sinhPhi/k;
        R(3,2) = k*sinhPhi;
    }

    // dispesions terms (only included if matrix is 6D)
    if(R.nrows()==6) {
        if(h && Kx)
            if(Kx>0) {
                double k = sqrt(Kx);
                double phi = k*length;
                double cosPhi = cos(phi);
                double sinPhi = sin(phi);

                R(0,5) = h*(1-cosPhi)/Kx;
                R(1,5) = h*sinPhi/k;
                R(4,5) = -h*h*(phi-sinPhi)/(Kx*k);
                R(4,0) = -h*sinPhi/k;
                R(4,1) = -h*(1-cosPhi)/Kx;
            }
            else {
                double k = sqrt(-Kx);
                double phi = k*length;
                double coshPhi = cosh(phi);
                double sinhPhi = sinh(phi);

                R(0,5) = h*(1-coshPhi)/Kx;
                R(1,5) = h*sinhPhi/k;
                R(4,5) = -h*h*(phi-sinhPhi)/(Kx*k);
                R(4,0) = -h*sinhPhi/k;
                R(4,1) = -h*(1-coshPhi)/Kx;
            }
    }
}

// solenoid R matrix including end field effects
// this matrix is symplectic

void SolR(double l, double ks, RealMatrix& rs)
{
    double theta=ks/2.*l;
    double c=cos(theta);
    double s=sin(theta);
    double c2=c*c;
    double sc=s*c;
    double s2=s*s;

    rs(0,0)=c2;
    rs(1,1)=c2;
    rs(2,2)=c2;
    rs(3,3)=c2;

    if(fequal(theta,0.)) {
        rs(0,1)=0;
        rs(2,3)=0;
    } else {
        rs(0,1)=2.*sc/ks;
        rs(2,3)=2.*sc/ks;
    }

    rs(0,2)=sc;
    rs(1,3)=sc;
    rs(2,0)=-sc;
    rs(3,1)=-sc;

    if(fequal(theta,0)) {
        rs(0,3)=0;
        rs(2,1)=0;
    }
    else {
        rs(0,3)= 2.*s2/ks;
        rs(2,1)=-2.*s2/ks;
    }

    rs(1,0)=-ks*sc/2;
    rs(3,2)=-ks*sc/2;

    rs(1,2)=-ks*s2/2;
    rs(3,0)= ks*s2/2;
}

}; // end of anonymous namespace


RealMatrix TransportMatrix::Drift (double length, RealMatrix& R)
{
    InitR(R);
    R(0,1)=length;
    if(R.nrows()>2) R(2,3)=length;
    return R;
}

RealMatrix& TransportMatrix::SectorBend (double length, double h, double Kx, RealMatrix& R)
{
    SBR(length,h,Kx,InitR(R));
    return R;
}

RealMatrix& TransportMatrix::SectorBendT (double l, double h, double K1, RealMatrix& T)
{
    // Currently only the chromatic part of the 6x6 matrix is supported
    InitR(T);

    if(l==0 || (K1==0 && h==0))
        return T;

    double h2 = h*h;
    double Kx = h2+K1;
    double Ky = -K1;
    double hk=2*h2+K1;

    double kx = Kx>0 ? sqrt(Kx) : sqrt(-Kx);
    double ky = Ky>0 ? sqrt(Ky) : sqrt(-Ky);

    // Horizontal plane
    if(Kx==0) {
        double l2=l*l;
        double l3=l*l2;
        T(0,0) = hk*l2/2.0;
        T(0,1) = hk*l3/6.0;
        T(1,0) = hk*l;
        T(1,1) = T(0,0);
    }
    else if(Kx>0) { // focusing
        double phix = l*kx;
        double cosphix = cos(phix);
        double sinphix = sin(phix);
        T(0,0)=T(1,1)=hk*l*sinphix/kx/2;
        T(0,1)=-hk*(phix*cosphix-sinphix)/(2*kx*kx*kx);
        T(1,0)=hk*(phix*cosphix+sinphix)/kx/2;
    }
    else { //defocusing
        double phix = l*kx;
        double coshphix = cos(phix);
        double sinhphix = sin(phix);
        T(0,0)=T(1,1)=hk*l*sinhphix/kx/2;
        T(0,1)=hk*(phix*coshphix-sinhphix)/(2*kx*kx*kx);
        T(1,0)=hk*(phix*coshphix+sinhphix)/kx/2;
    }

    // Vertical plane (as for QuadrupoleT)
    if(Ky>0) {// focusing
        double phiy = l*ky;
        double cosphiy = cos(phiy);
        double sinphiy = sin(phiy);
        T(2,2) = T(3,3) = phiy*sinphiy/2;
        T(2,3) = (sinphiy/ky-l*cosphiy)/2;
        T(3,2) = ky*(phiy*cosphiy+sinphiy)/2;
    }
    else if(Ky<0) {// defocusing
        double phiy = l*ky;
        double coshphiy = cosh(phiy);
        double sinhphiy = sinh(phiy);
        T(2,2) = T(3,3) = -phiy*sinhphiy/2;
        T(2,3) = (sinhphiy/ky-l*coshphiy)/2;
        T(3,2) = -ky*(phiy*coshphiy+sinhphiy)/2;
    }

    // Dispersion/time terms
    if(Kx==0) {
        double phi = h*l;
        double phi2=phi*phi;
        T(0,5) = phi*l*(phi2-12.0)/24.0;
        T(1,5) = phi*(phi2/6.0-1.0);
        T(4,5) = -phi2*phi*(phi2-40.0)/120.0;
    } else if(Kx>0) {
        double phi = l*kx;
        double sinphi = sin(phi);
        double cosphi = cos(phi);
        double Kx2=Kx*Kx;
        T(0,5) = -h*(2*h2*(cosphi-1)+hk*phi*sinphi)/(2*Kx2);
        T(1,5) = -h*(hk*phi*cosphi+K1*sinphi)/(2*pow(kx,3));
        T(4,5) =  h2*(2*K1*phi-hk*phi*cosphi+(2*h2-K1)*sinphi)/(2*pow(kx,5));
    } else {
        double phi = l*kx;
        double sinhphi = sinh(phi);
        double coshphi = cosh(phi);
        double Kx2=Kx*Kx;
        T(0,5) = h*(2*h2*(1-coshphi)+hk*phi*sinhphi)/(2*Kx2);
        T(1,5) = h*(hk*phi*coshphi+K1*sinhphi)/(2*pow(kx,3));
        T(4,5) = h2*(2*K1*phi-hk*phi*coshphi+(2*h2-K1)*sinhphi)/(2*pow(kx,5));
    }

    return T;
}

RealMatrix& TransportMatrix::QuadrupoleR (double length, double Kx, RealMatrix& R)
{
    SBR(length,0,Kx,InitR(R));
    return R;
}

RealMatrix& TransportMatrix::QuadrupoleT (double length, double Kx, RealMatrix& T)
{
    InitR(T);

    if(length==0 || Kx==0)
        return T;

    double k = Kx>0 ? sqrt(Kx) : sqrt(-Kx);
    double phi = length*k;
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double coshphi = cosh(phi);
    double sinhphi = sinh(phi);


    T(0,0) = T(1,1) = phi*sinphi/2;
    T(0,1) = (sinphi/k-length*cosphi)/2;
    T(1,0) = k*(phi*cosphi+sinphi)/2;

    T(2,2) = T(3,3) = -phi*sinhphi/2;
    T(2,3) = (sinhphi/k-length*coshphi)/2;
    T(3,2) = -k*(phi*coshphi+sinhphi)/2;

    if(Kx<0) { // defocusing
        swap(T(0,0),T(2,2));
        swap(T(0,1),T(2,3));
        swap(T(1,0),T(3,2));
        swap(T(1,1),T(3,3));
    }

    return T;
}

RealMatrix& TransportMatrix::Srot (double cosphi, double sinphi, RealMatrix& R)
{
    InitR(R);
    SROT(cosphi,sinphi,R);
    return R;
}

RealMatrix& TransportMatrix::SrotR4 (double cosphi, double sinphi, RealMatrix& R)
{
    InitR(R);
    SROT(cosphi,sinphi,R);
    return R;
}

RealMatrix& TransportMatrix::PoleFaceRot (double h, double theta, double fint, double hgap, RealMatrix& R)
{
    InitR(R);
    const double sinTheta = sin(theta);
    const double phi = 2.0*fint*hgap*h*(1+sinTheta*sinTheta)/cos(theta);
    R(1,0) = h*tan(theta);
    R(3,2) =-h*tan(theta-phi);
    return R;
}


RealMatrix& TransportMatrix::Solenoid (double length, double K0, double K1, bool entrFringe, bool exitFringe, RealMatrix& R)
{
    InitR(R);
    SolR(length,K0,R);
    return R;

    /**********
    double K = K0/2.0;
    double phi = K*length;
    double C = cos(phi);
    double S = sin(phi);

    R(0,0)=R(1,1)=R(2,2)=R(3,3)=C*C;

    R(0,1)= S*C/K;
    R(0,2)= S*C;
    R(0,3)= S*S/K;

    R(1,0)=-K*S*C;
    R(1,2)=-K*S*S;
    R(1,3)= S*C;

    R(2,0)=-S*C;
    R(2,1)=-S*S/K;
    R(2,3)= S*C/K;

    R(3,0)= K*S*S;
    R(3,1)=-S*C;
    R(3,2)=-K*S*C;

    return R;
    ***/
}

RealMatrix& TransportMatrix::TWRFCavity (double length, double g, double f, double phi,
        double E0, bool inc_end_fields, RealMatrix& R)
{
    using namespace PhysicalConstants;

    InitR(R);

    const double dE=length*g;
    const double cosphi = cos(phi);
    const double dEcosPhi = dE*cosphi;
    const double E1 = E0+dEcosPhi;
    const double Er = E1/E0;
    const double logEr = log(Er);

    // Condition for zero acceleration modified by A.Wolski, 5 November 2003
    // Loss of numerical accuracy can lead to dEcosPhi of order 1e-20
    // in which case the transport matrix can have large erroneous terms.
    if(fabs(dEcosPhi)<1.0e-16) {
        R(0,1)=length;
    }
    else {
        if(!inc_end_fields) {
            R(0,0) = 1-0.5*logEr;
            R(0,1) = E0*length*logEr/dEcosPhi;
            R(1,0) = -dEcosPhi*logEr/(4*E1*length);
            R(1,1) = (1+0.5*logEr)/Er;
        }
        else {
            R(0,1) = length*E0*logEr/dEcosPhi;
            R(1,1) = Er;
        }
    }

    if(R.nrows()>2) {
        R(2,2)=R(0,0);
        R(2,3)=R(0,1);
        R(3,2)=R(1,0);
        R(3,3)=R(1,1);
    }

    if(R.nrows()==6) {
        const double k = twoPi*f/SpeedOfLight;
        R(5,4)=k*dE*sin(phi)/E1;
        R(5,5)=1/Er;
    }
    return R;
}

RealMatrix& TransportMatrix::SWRFCavity (int ncells, double g, double f, double phi, double E0, RealMatrix& R)
{
    using namespace PhysicalConstants;

    InitR(R);

    double root2 = sqrt(2.0);
    double root8 = sqrt(8.0);

    double lambda = SpeedOfLight/f;
    double length = ncells*lambda/2;

    double cosPhi = cos(phi);
    double E1=E0+g*length*cosPhi;

    double alpha = log(E1/E0)/cosPhi/root8;
    double cosAlpha = cos(alpha);
    double sinAlpha = sin(alpha);


    R(0,0) = cosAlpha - root2*cosPhi*sinAlpha;
    R(0,1) = root8*E0*sinAlpha/g;
    R(1,0) = -g*(2+cos(2*phi))*sinAlpha/E1/root8;
    R(1,1) = E0*(cosAlpha+root2*cosPhi*sinAlpha)/E1;

    if(R.nrows()==4) {
        R(2,2)=R(0,0);
        R(2,3)=R(0,1);
        R(3,2)=R(1,0);
        R(3,3)=R(1,1);
    }

    return R;
}

RealMatrix& TransportMatrix::Arb3DXForm (const RealMatrix& R3d, const RealVector& X, RealMatrix& R, RealVector& dX)
{
    // linear approximation for a general 3-D Euclidean transformation
    // for TRANSPORT canonical variables.
    // Given a general 3-D rotation matrix r (3x3) and a translation u, then we can
    // calculate a linear transformation (R matrix) and a 6-D displacement (X).
    // Note that R is in general not symplectic, and large x- or y-axis rotations
    // will result in errors.

    // Note that the expressions where derived
    // (and generated) using Mathematica.

    const double r11 = R3d(0,0);
    const double r12 = R3d(0,1);
    const double r13 = R3d(0,2);
    const double r21 = R3d(1,0);
    const double r22 = R3d(1,1);
    const double r23 = R3d(1,2);
    const double r31 = R3d(2,0);
    const double r32 = R3d(2,1);
    const double r33 = R3d(2,2);
    const double x = X[0];
    const double y = X[1];
    const double z = X[2];

    const double A1 = 1/r33;
    const double A0 = A1*A1;
    const double A2 = -(r31*x);
    const double A3 = -(r32*y);
    const double A4 = -(r33*z);
    const double A5 = A2 + A3 + A4;
    const double A6 = A1*A5;

    R.copy(IdentityMatrix(6));

    R(0,0) = r11 - A1*r13*r31;
    R(0,1) = -(A1*A5*r11) + A0*A5*r13*r31;
    R(0,2) = r12 - A1*r13*r32;
    R(0,3) = -(A1*A5*r12) + A0*A5*r13*r32;
    R(1,1) = A1*r11 - A0*r13*r31;
    R(1,3) = A1*r12 - A0*r13*r32;
    R(2,0) = r21 - A1*r23*r31;
    R(2,1) = -(A1*A5*r21) + A0*A5*r23*r31;
    R(2,2) = r22 - A1*r23*r32;
    R(2,3) = -(A1*A5*r22) + A0*A5*r23*r32;
    R(3,1) = A1*r21 - A0*r23*r31;
    R(3,3) = A1*r22 - A0*r23*r32;
    R(4,0) = r31;
    R(4,2) = r32;
    R(4,4) = r33;

    dX.redim(6);

    dX[0]= -r11*x-r12*y-r13*z-r13*A6;
    dX[1]=  A1*r13;
    dX[2]= -r21*x-r22*y-r23*z-r23*A6;
    dX[3]=  A1*r23;
    dX[4]=  A5;
    dX[5]=  0;

    return R;
}

