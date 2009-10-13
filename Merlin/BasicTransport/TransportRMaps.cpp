/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/10/24 09:55:39 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/MatrixPrinter.h"
#include "NumericalUtils/utils.h"
#include "BasicTransport/TransportRMap.h"

namespace TransportRMap {

inline double LogOnePlusX(double x) { return x*(1-0.5*x); }

void SROT(double cosphi, double sinphi, RMap& R)
{
    R.AddTerm(1,1,cosphi);
    R.AddTerm(2,2,cosphi);
    R.AddTerm(3,3,cosphi);
    R.AddTerm(4,4,cosphi);
    R.AddTerm(1,3,-sinphi);
    R.AddTerm(2,4,-sinphi);
    R.AddTerm(3,1,sinphi);
    R.AddTerm(4,2,sinphi);
}

void QR2(double length, double K1, int offset, RMap& R)
{
    if(K1==0) {
        R.AddTerm(offset+1,offset+1,1.0);
        R.AddTerm(offset+1,offset+2,length);
        R.AddTerm(offset+2,offset+2,1.0);
    }
    else if(K1>0) {
        //  focusing magnet
        double k = sqrt(K1);
        double phi = k*length;
        double sinphi = sin(phi);
        double cosphi = cos(phi);
        R.AddTerm(offset+1,offset+1,cosphi);
        R.AddTerm(offset+2,offset+2,cosphi);
        R.AddTerm(offset+1,offset+2,sinphi/k);
        R.AddTerm(offset+2,offset+1,-k*sinphi);
    }
    else {
        //  defocusing magnet
        double k = sqrt(-K1);
        double phi = k*length;
        double sinhphi = sinh(phi);
        double coshphi = cosh(phi);
        R.AddTerm(offset+1,offset+1,coshphi);
        R.AddTerm(offset+2,offset+2,coshphi);
        R.AddTerm(offset+1,offset+2,sinhphi/k);
        R.AddTerm(offset+2,offset+1,k*sinhphi);
    }
}

void SBR(double length, double h, double Kx, RMap& R)
{
    // check for zero length element
    if(length==0) {
        MakeIdentity(R);
        return;
    }

    double Ky = -Kx;
    Kx += h*h;

    // horizontal plane
    QR2(length,Kx,0,R);

    // vectircal plane
    QR2(length,Ky,2,R);

    if(h && Kx) {
        if(Kx>0) {
            double k = sqrt(Kx);
            double phi = k*length;
            double cosPhi = cos(phi);
            double sinPhi = sin(phi);

            R.AddTerm(1,6) = h*(1-cosPhi)/Kx;
            R.AddTerm(2,6) = h*sinPhi/k;
            R.AddTerm(5,6) = -h*h*(phi-sinPhi)/(Kx*k);
            R.AddTerm(5,1) = -h*sinPhi/k;
            R.AddTerm(5,2) = -h*(1-cosPhi)/Kx;
        }
        else {
            double k = sqrt(-Kx);
            double phi = k*length;
            double coshPhi = cosh(phi);
            double sinhPhi = sinh(phi);

            R.AddTerm(1,6) = h*(1-coshPhi)/Kx;
            R.AddTerm(2,6) = h*sinhPhi/k;
            R.AddTerm(5,6) = -h*h*(phi-sinhPhi)/(Kx*k);
            R.AddTerm(5,1) = -h*sinhPhi/k;
            R.AddTerm(5,2) = -h*(1-coshPhi)/Kx;
        }
    }

    R.AddTerm(5,5)=1.0;
    R.AddTerm(6,6)=1.0;
}

// solenoid R matrix including end field effects
// this matrix is symplectic

void SolR(double l, double ks, RMap& rs)
{
    double theta=ks/2.*l;
    double c=cos(theta);
    double s=sin(theta);
    double c2=c*c;
    double sc=s*c;
    double s2=s*s;

    rs.AddTerm(1,1)=c2;
    rs.AddTerm(2,2)=c2;
    rs.AddTerm(3,3)=c2;
    rs.AddTerm(4,4)=c2;

    if(!fequal(theta,0.)) {
        rs.AddTerm(1,2)=2.*sc/ks;
        rs.AddTerm(3,4)=2.*sc/ks;
    }

    rs.AddTerm(1,3)=sc;
    rs.AddTerm(2,4)=sc;
    rs.AddTerm(3,1)=-sc;
    rs.AddTerm(4,2)=-sc;

    if(!fequal(theta,0)) {
        rs.AddTerm(1,4)= 2.*s2/ks;
        rs.AddTerm(4,2)=-2.*s2/ks;
    }

    rs.AddTerm(2,1)=-ks*sc/2;
    rs.AddTerm(4,3)=-ks*sc/2;

    rs.AddTerm(2,3)=-ks*s2/2;
    rs.AddTerm(4,1)= ks*s2/2;
}



// Class Utility TransportRMap


void Drift (double length, RMap& R)
{
    MakeIdentity(R);
    if(length!=0) {
        R.AddTerm(1,2,length);
        R.AddTerm(3,4,length);
    }
}

void SectorBend(double length, double h, double Kx, RMap& R)
{
    SBR(length,h,Kx,R);
}

void Quadrupole(double length, double Kx, RMap& R)
{
    QR2(length,Kx,0,R);
    QR2(length,-Kx,2,R);
    //		R.AddTerm(5,5,1.0);
    //		R.AddTerm(6,6,1.0);
}

void Srot(double cosphi, double sinphi, RMap& R)
{
    SROT(cosphi,sinphi,R);
}


void PoleFaceRot(double h, double theta, double fint, double hgap, RMap& R)
{

    MakeIdentity(R);
    const double sinTheta = sin(theta);
    const double phi = 2.0*fint*hgap*h*(1+sinTheta*sinTheta)/cos(theta);
    R.AddTerm(2,1) = h*tan(theta);
    R.AddTerm(4,3) =-h*tan(theta-phi);
}


void Solenoid (double length, double K0, double K1, bool entrFringe, bool exitFringe, RMap& R)
{
    SolR(length,K0,R);
}

void TWRFCavity (double length, double g, double f, double phi,
                 double E0, bool entr_field, bool exit_field, RMap& R)
{

    using namespace PhysicalConstants;

    const double dE=length*g;
    const double cosphi = cos(phi);
    const double dEcosPhi = dE*cosphi;
    const double E1 = E0+dEcosPhi;
    const double Er = E0/E1;
    const double logEr = LogOnePlusX(dEcosPhi/E0);

    if(dEcosPhi==0) {
        MakeIdentity(R);
        R.AddTerm(1,2)=R.AddTerm(3,4)=length;
    }
    else {
        if(entr_field && !exit_field) {
            R.AddTerm(1,1) = R.AddTerm(3,3) = 1-0.5*logEr/cosphi;
            R.AddTerm(1,2) = R.AddTerm(3,4) = E0*logEr/g/cosphi;
            R.AddTerm(2,1) = R.AddTerm(4,3) =-0.5*g/E1;
            R.AddTerm(2,2) = R.AddTerm(4,4) = Er;
        }
        else if(!entr_field && exit_field) {
            R.AddTerm(1,1) = R.AddTerm(3,3) = 1;
            R.AddTerm(1,2) = R.AddTerm(3,4) = E0*logEr/g/cosphi;
            R.AddTerm(2,1) = R.AddTerm(4,3) = 0.5*g*cosphi/E1;
            R.AddTerm(2,2) = R.AddTerm(4,4) = Er*(1+0.5*logEr/cosphi);
        }
        else if(entr_field && exit_field) {
            R.AddTerm(1,1) = R.AddTerm(3,3) = 1-0.5*logEr/cosphi;
            R.AddTerm(1,2) = R.AddTerm(3,4) = E0*logEr/g/cosphi;
            R.AddTerm(2,1) = R.AddTerm(4,3) = -0.25*g*logEr/E1/cosphi;
            R.AddTerm(2,2) = R.AddTerm(4,4) = Er*(1+0.5*logEr/cosphi);
        }
        else { //no fringe fields
            R.AddTerm(1,1) = R.AddTerm(3,3) = 1.0;
            R.AddTerm(1,2) = R.AddTerm(3,4) = E0*logEr/g/cosphi;
            R.AddTerm(2,2) = R.AddTerm(4,4) = Er;
        }
    }
}

void SWRFCavity (int ncells, double g, double f, double phi, double E0, RMap& R)
{

    using namespace PhysicalConstants;

    double root2 = sqrt(2.0);
    double root8 = sqrt(8.0);

    double lambda = SpeedOfLight/f;
    double length = ncells*lambda/2;

    double cosPhi = cos(phi);
    double E1=E0+g*length*cosPhi;

    double alpha = log(E1/E0)/cosPhi/root8;
    double cosAlpha = cos(alpha);
    double sinAlpha = sin(alpha);

    R.AddTerm(1,1) = R.AddTerm(3,3) = cosAlpha - root2*cosPhi*sinAlpha;
    R.AddTerm(1,2) = R.AddTerm(3,4) = root8*E0*sinAlpha/g;
    R.AddTerm(2,1) = R.AddTerm(4,3) = -g*(2+cos(2*phi))*sinAlpha/E1/root8;
    R.AddTerm(2,2) = R.AddTerm(4,4) = E0*(cosAlpha+root2*cosPhi*sinAlpha)/E1;
}

// functions returning R2Map objects
void Drift (double length, R2Map& R)
{
    R.r11=R.r22=1.0;
    R.r12=length;
    R.r21=0;
}

void Quadrupole (double length, double K1, R2Map& R)
{
    if(K1==0)
        Drift(length,R);
    else if(K1>0) {
        //  focusing magnet
        double k = sqrt(K1);
        double phi = k*length;
        double sinphi = sin(phi);
        double cosphi = cos(phi);
        R.r11=R.r22=cosphi;
        R.r12=sinphi/k;
        R.r21=-k*sinphi;
    }
    else {
        //  defocusing magnet
        double k = sqrt(-K1);
        double phi = k*length;
        double sinhphi = sinh(phi);
        double coshphi = cosh(phi);
        R.r11=R.r22=coshphi;
        R.r12=sinhphi/k;
        R.r21=k*sinhphi;
    }
}

void TWRFCavity (double length, double g, double f, double phi,
                 double E0, bool entr_field, bool exit_field, R2Map& R)
{
    using namespace PhysicalConstants;

    const double dE=length*g;
    const double cosphi = cos(phi);
    const double dEcosPhi = dE*cosphi;
    const double E1 = E0+dEcosPhi;
    const double Er = E0/E1;
    const double logEr = LogOnePlusX(dEcosPhi/E0);
    /****
    		if(dEcosPhi==0)
    		Drift(length,R);
    	else {
    		if(entr_field && !exit_field) {
    			R.r11= 1-0.5*logEr/cosphi;
    			R.r12= E0*logEr/g/cosphi;
    			R.r21=-0.5*g/E1;
    			R.r22= Er;
    		}
    		else if(!entr_field && exit_field) {
    			R.r11= 1;
    			R.r12= E0*logEr/g/cosphi;
    			R.r21= 0.5*g/E1;
    			R.r22= Er*(1+0.5*logEr/cosphi);
    		}
    		else if(entr_field && exit_field) {
    			R.r11= 1-0.5*logEr/cosphi;
    			R.r12= E0*logEr/g/cosphi;
    			R.r21= -0.25*g*logEr/E1/cosphi;
    			R.r22= Er*(1+0.5*logEr/cosphi);
    		}
    		else { // no fringe fields
    			R.r11= 1.0;
    			R.r12= E0*logEr/g/cosphi;
    			R.r21= 0;
    			R.r22= Er;
    		}
    	}
    ****/
    if(dEcosPhi==0)
        Drift(length,R);
    else {
        if(entr_field && !exit_field) {
            R.r11= 1-0.5*logEr;
            R.r12= E0*logEr/g/cosphi;
            R.r21=-0.5*g*cosphi/E1;
            R.r22= Er;
        }
        else if(!entr_field && exit_field) {
            R.r11= 1;
            R.r12= E0*logEr/g/cosphi;
            R.r21= 0.5*g*cosphi/E1;
            R.r22= Er*(1+0.5*logEr);
        }
        else if(entr_field && exit_field) {
            R.r11= 1-0.5*logEr;
            R.r12= E0*logEr/g/cosphi;
            R.r21= -0.25*g*logEr*cosphi/E1;
            R.r22= Er*(1+0.5*logEr);
        }
        else { // no fringe fields
            R.r11= 1.0;
            R.r12= E0*logEr/g/cosphi;
            R.r21= 0;
            R.r22= Er;
        }
    }
}

void SWRFCavity (int ncells, double g, double f, double phi, double E0, R2Map& R)
{
    using namespace PhysicalConstants;

    double root2 = sqrt(2.0);
    double root8 = sqrt(8.0);

    double lambda = SpeedOfLight/f;
    double length = ncells*lambda/2;

    double cosPhi = cos(phi);
    double E1=E0+g*length*cosPhi;

    double alpha = log(E1/E0)/cosPhi/root8;
    double cosAlpha = cos(alpha);
    double sinAlpha = sin(alpha);

    R.r11 = cosAlpha - root2*cosPhi*sinAlpha;
    R.r12 = root8*E0*sinAlpha/g;
    R.r21 = -g*(2+cos(2*phi))*sinAlpha/E1/root8;
    R.r22 = E0*(cosAlpha+root2*cosPhi*sinAlpha)/E1;
}

}; // end namespace TransportRMap
