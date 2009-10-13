/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/26 19:30:15 $
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include "BasicTransport/NormalTransform.h"
#include "BasicTransport/RMap.h"
#include <fstream>

double ProjectedEmittance(const PSmoments& s, PScoord x1, PScoord x2)
{
    return sqrt(s.var(x1)*s.var(x2)-pow(s(x1,x2),2));
}

RealMatrix CouplingMatrix(double a, double b, double c, double d)
{
    RealMatrix C=IdentityMatrix(6);

    const double detB2 = b*c-a*d;
    double detB;
    double h;
    double zeta;

    if(detB2==0)
        return C;

    if(detB2>0) {
        detB = sqrt(detB2);
        h=cos(detB);
        zeta=sin(detB)/detB;
    }
    else{
        detB = sqrt(-detB2);
        h=cosh(detB);
        zeta=sinh(detB)/detB;
    }

    a*=zeta;
    b*=zeta;
    c*=zeta;
    d*=zeta;

    C(0,0)=h;
    C(1,1)=h;
    C(2,2)=h;
    C(3,3)=h;

    C(0,2)=a;
    C(0,3)=b;
    C(1,2)=-c;
    C(1,3)=-d;

    C(2,0)=d;
    C(2,1)=b;
    C(3,0)=-c;
    C(3,1)=-a;

    return C;
}

RealMatrix AlphaMatrix(double ax, double ay)
{
    RealMatrix A = IdentityMatrix(6);
    A(1,0)=-ax;
    A(3,2)=-ay;
    return A;
}

RealMatrix BetaMatrix(double bx, double by)
{
    RealMatrix B = IdentityMatrix(6);
    const double sbx=sqrt(bx);
    const double sby=sqrt(by);
    B(0,0)=sbx;
    B(1,1)=1/sbx;
    B(2,2)=sby;
    B(3,3)=1/sby;
    return B;
}

RealMatrix DispersionMatrix(double Dx, double Dxp, double Dy, double Dyp)
{
    RealMatrix D=IdentityMatrix(6);
    D(0,5)=Dx;
    D(1,5)=Dxp;
    D(2,5)=Dy;
    D(3,5)=Dyp;
    return D;
}

RealMatrix InverseBetaTransform(double bx, double by, double ax, double ay)
{
	RealMatrix R=IdentityMatrix(6);
	double sqrtBX = sqrt(bx);
	double sqrtBY = sqrt(by);

	R(0,0)=1/sqrtBX;
	R(1,0)=ax/sqrtBX;
	R(1,1)=sqrtBX;

	R(2,2)=1/sqrtBY;
	R(3,2)=ay/sqrtBY;
	R(3,3)=sqrtBY;

	return R;
}

RealMatrix NormalTransform(const BeamData& t)
{
    //  assert(t.ok());
    RealMatrix D=DispersionMatrix(t.Dx,t.Dxp,t.Dy,t.Dyp);
    RealMatrix B=BetaMatrix(t.beta_x,t.beta_y);
    RealMatrix A=AlphaMatrix(t.alpha_x,t.alpha_y);
    RealMatrix C=CouplingMatrix(t.c_xy,t.c_xyp,t.c_xpy,t.c_xpyp);

    return D*B*A*C;
}


PSmoments& BeamDataToSigmaMtrx(const BeamData& t, PSmoments& S)
{
    RMap R(NormalTransform(t));
    S.zero();
    S(0,0)=t.emit_x;
    S(1,1)=t.emit_x;
    S(2,2)=t.emit_y;
    S(3,3)=t.emit_y;
    S(4,4)=pow(t.sig_z,2);
    S(5,5)=pow(t.sig_dp,2);
    R.Apply(S);

    S[0]=t.x0;
    S[1]=t.xp0;
    S[2]=t.y0;
    S[3]=t.yp0;

    return S;
}

RealMatrix DecoupleSigma(SigmaMatrix& S)
{
	// This routine removes the x-y coupling from S0 and returns the 
	// corresponding C matrix.
	RealVector c(1.0,4); // (a,b,c,d)
	RealMatrix R = IdentityMatrix(6);

	do{
		double d = S(2,2)*S(3,3)-S(2,3)*S(2,3)-(S(0,0)*S(1,1)-S(0,1)*S(0,1));
		c(0)=(-S(0,1)*S(0,3)+S(0,0)*S(1,3)-S(0,3)*S(2,3)+S(0,2)*S(3,3))/d;
		c(1)=(S(0,1)*S(0,2)-S(0,0)*S(1,2)+S(0,3)*S(2,2)-S(0,2)*S(2,3))/d;
		c(2)=(S(0,3)*S(1,1)-S(0,1)*S(1,3)+S(1,3)*S(2,3)-S(1,2)*S(3,3))/d;
		c(3)=(-S(0,2)*S(1,1)+S(0,1)*S(1,2)-S(1,3)*S(2,2)+S(1,2)*S(2,3))/d;
		RealMatrix R1 = CouplingMatrix(-c(0),-c(1),-c(2),-c(3));
		RMap(R1).Apply(S);
		R=R1*R;
	} while(c*c>=1.0e-10);

	return R;
}


BeamData& SigmaMatrixToBeamData(const PSmoments& S0, BeamData& t)
{
    PSmoments S=S0;
    t=BeamData();

    // First check if we have dispersion:
    double d2 = S.var(ps_DP);
    if(d2!=0) {
        t.Dx = S(ps_X,ps_DP)/d2;
        t.Dxp = S(ps_XP,ps_DP)/d2;
        t.Dy = S(ps_Y,ps_DP)/d2;
        t.Dyp = S(ps_YP,ps_DP)/d2;

        // Remove energy correlations
        RMap D(DispersionMatrix(-t.Dx,-t.Dxp,-t.Dy,-t.Dyp));
        D.Apply(S);
    }

	// Remove the cross-plane coupling
	// and calculate the coupling parameters
	RealMatrix C=DecoupleSigma(S);
	double z;
	if(C(0,0)>1.0) {
		double x=sqrt(C(0,0)*C(0,0)-1);
		z=x!=0 ? log(C(0,0)+x)/x : 0;
	}
	else {
		double sinPhi = sqrt(1-C(0,0));
		double phi = atan2(sinPhi,C(0,0));
		z= sinPhi!=0 ? phi/sinPhi : 0;
	}
	t.c_xy  = C(0,2)*z;
	t.c_xyp = C(0,3)*z;
	t.c_xpy = C(1,2)*z;
	t.c_xpyp= C(1,3)*z;

    t.emit_x = ProjectedEmittance(S,ps_X,ps_XP);
    t.emit_y = ProjectedEmittance(S,ps_Y,ps_YP);

	// Now calculate beta, alpha, 
    t.beta_x = S.var(ps_X)/t.emit_x;
    t.beta_y = S.var(ps_Y)/t.emit_y;
    t.alpha_x = -S(ps_X,ps_XP)/t.emit_x;
    t.alpha_y = -S(ps_Y,ps_YP)/t.emit_y;

	// Remove the correlation
	RealMatrix B=InverseBetaTransform(t.beta_x,t.beta_y,t.alpha_x,t.alpha_y);
 	RMap(B).Apply(S);

    t.sig_z  = S0.std(ps_CT);
    t.sig_dp = S0.std(ps_DP);

    // Centroid
    t.x0  = S0[0];
    t.xp0 = S0[1];
    t.y0  = S0[2];
    t.yp0 = S0[3];
    t.ct0 = S0[4];
    t.p0  = S0[5]; // not actually the energy, but dp/p!

    return t;
}

pair<double,double> NormalModeEmittance(const PSmoments& S)
{
	BeamData t;
	SigmaMatrixToBeamData(S,t);
	return pair<double,double>(t.emit_x,t.emit_y);
}
