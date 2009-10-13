/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "BeamDynamics/ParticleTracking/Integrators/LCAVintegrator.h"
#include "BasicTransport/BasicTransportMaps.h"
#include "NumericalUtils/PhysicalConstants.h"

using namespace PhysicalConstants;

#define CHK_ZERO(s) if(s==0) return;

namespace ParticleTracking {

struct EntranceFieldMap {

    double k;
    double Ez;
    double phi0;

    EntranceFieldMap(double g, double k1, double phi, double p)
            : Ez(g/p),k(k1),phi0(phi) {}

    void Apply(PSvector& x) const {
        double a = -0.5*Ez*cos(phi0-k*x.ct())/(1+x.dp());
        x.xp() += a*x.x();
        x.yp() += a*x.y();
    }
};

struct DriftMap {
    double s;
    DriftMap(double l) : s(l) {}

    void Apply(PSvector& x) const {
        x.x()+=x.xp()*s;
        x.y()+=x.yp()*s;
    }
};

struct LCAVMap {

    double k;
    double Ez;
    double phi0;
    double L;
    double E0;
    double E1;
    mutable double Esum;
    mutable size_t np;

    LCAVMap(double g, double ds, double k1, double phi, double p0)
            : k(k1),Ez(g*ds),L(ds),phi0(phi),E0(p0),E1(p0+g*ds*cos(phi)),Esum(0),np(0)
    {}

    void Apply(PSvector& x) const {

        double cosphi = cos(phi0-k*x.ct());
        double dE = Ez*cosphi;
        double Ein = E0*(1+x.dp());
        double Eout = Ein+dE;
        double fact = Ein/Eout;
        //			double r12 = L*Ein*log(1+dE/Ein)/dE;
        double r12 = L*(1-0.5*dE/Ein); // 2nd-order expansion of log(1+x)

        x.x()+=r12*x.xp();
        x.xp()*=fact;
        x.y()+=r12*x.yp();
        x.yp()*=fact;

        x.dp() = Eout/E1-1.0;

        Esum += Eout;
        np++;
    }

    double Eav() const { return E1; }
};

// Functor ApplyRFdp (used for no change in reference momentum)
struct ApplyRFMap {

    double Vn,k,phi0,d0,ds;

    ApplyRFMap(double Vnorm, double kval, double phase, double len)
            : Vn(Vnorm),k(kval),phi0(phase),ds(len) {};

    void Apply(PSvector& p) const {
        double ddp = Vn*cos(phi0-k*p.ct());
        p.dp() += ddp;
        p.x()  += ds*p.xp();
        p.y()  += ds*p.yp();
    }
};

void LCAVIntegrator::TrackEntrance ()
{
    ApplyEndField(1.0);
}

void LCAVIntegrator::TrackExit ()
{
    ApplyEndField(-1.0);
}

void LCAVIntegrator::ApplyEndField(double gsgn)
{
    const RFAcceleratingField& field = (*currentComponent).GetField();
    double g   = field.GetAmplitude();
    double k   = field.GetK();
    double phi = field.GetPhase();
    double E0  = currentBunch->GetReferenceMomentum();
    if(g!=0)
        ApplyMap(EntranceFieldMap(gsgn*g,k,phi,E0),currentBunch->GetParticles());
}

void LCAVIntegrator::TrackStep (double ds)
{
    CHK_ZERO(ds);

    const RFAcceleratingField& field = (*currentComponent).GetField();

    double g   = field.GetAmplitude();
    double k   = field.GetK();
    double phi = field.GetPhase();
    bool full_acceleration = field.FullAcceleration();
    double E0  = currentBunch->GetReferenceMomentum();

    if(g==0) {
        ApplyMap(DriftMap(ds),currentBunch->GetParticles());
        return;
    }


    if(full_acceleration) {
        // structure map
        LCAVMap lcmap(g,ds,k,phi,E0);
        ApplyMap(lcmap,currentBunch->GetParticles());
        currentBunch->SetReferenceMomentum(lcmap.Eav());
    } else {
        ApplyMap(ApplyRFMap(g*ds/E0,k,phi,ds),currentBunch->GetParticles());
    };

    return;
}

}; // end namespace ParticleTracking


