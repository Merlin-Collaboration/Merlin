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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "BeamDynamics/ParticleTracking/Integrators/TransRFIntegrator.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "BasicTransport/RMap.h"

using namespace PhysicalConstants;

#define CHK_ZERO(s) if(s==0) return;

namespace ParticleTracking {

struct TRFMap {
    
    double cosTheta,sinTheta,gRed;
    double k, phi0, len;

    TRFMap(double g, double ds, double k1, double phi, double theta, double p0)
        :cosTheta(cos(theta)),sinTheta(sin(theta)),gRed(g/p0),k(k1),phi0(phi),len(ds)
    {}

    void Apply(PSvector& x) const {

        double dA = gRed*len*cos(phi0-k*x.ct())/(1+x.dp());
        double cx1  = cosTheta*dA;
        double cy1  = sinTheta*dA;
        double cx2  = 0.5*cx1*len;
        double cy2  = 0.5*cy1*len;

        x.xp()+=cx1;
        x.x() +=cx2+len*x.xp();
        x.yp()+=cy1;
        x.y() +=cy2+len*x.yp();
    }
};


void TransRFIntegrator::TrackEntrance ()
{
    // TO DO
}

void TransRFIntegrator::TrackExit ()
{
    // TO DO
}

void TransRFIntegrator::TrackStep (double ds)
{
    CHK_ZERO(ds);

    const TransverseRFfield& field = static_cast<const TransverseRFfield&>(currentComponent->GetField());

    double g   = field.GetAmplitude();
    double k   = field.GetK();
    double phi = field.GetPhase();
    double E0  = currentBunch->GetReferenceMomentum();
    double theta = field.GetFieldOrientation();

    if(g==0)
        ApplyMap(ApplyDriftWithPathLength(ds),currentBunch->GetParticles());
    else    
        ApplyMap(TRFMap(g,ds,k,phi,theta,E0),currentBunch->GetParticles());

    return;
}

}; // end namespace ParticleTracking


