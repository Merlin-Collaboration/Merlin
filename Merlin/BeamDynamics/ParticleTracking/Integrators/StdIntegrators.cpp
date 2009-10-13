/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/26 20:12:15 $
// $Revision: 1.9 $
// 
/////////////////////////////////////////////////////////////////////////

#include <limits>
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/MatrixPrinter.h"
#include "BasicTransport/TransportMatrix.h"
#include "BasicTransport/MatrixMaps.h"
#include "BeamDynamics/ParticleTracking/Integrators/StdIntegrators.h"
#include "BeamDynamics/ParticleTracking/Integrators/LCAVintegrator.h"
#include "BeamDynamics/ParticleTracking/Integrators/TransRFIntegrator.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

#define CHK_ZERO(s) if(s==0) return

namespace {

// Functor which applies a drift
struct psdrift {
    double z;
    psdrift(double len):z(len){}
    void operator()(PSvector& p) {
        p.x() += p.xp()*z;
        p.y() += p.yp()*z;
    }
};

inline void ApplyDrift(PSvectorArray& psv,double z)
{
    if(z!=0)
        for_each(psv.begin(),psv.end(),psdrift(z));
}

struct MultipoleKick {
    const MultipoleField& field;
    double scale;

    MultipoleKick(const MultipoleField& f, double len, double P0, double q)
            : field(f),scale(q*len*eV*SpeedOfLight/P0) {}

    void operator()(PSvector& v) {
        double x=v.x();
        double y=v.y();
        double dp=v.dp();
        Complex F = scale*field.GetField2D(x,y)/(1+dp);
        v.xp() += -F.real();
        v.yp() +=  F.imag();
    };
};


// Functor ApplyRFdp (used for full acceleration)
struct ApplyRFdp {

    double Vn,k,phi0,cosPhi0,d0;
    const RMtrx& R;
    bool fullacc;

    ApplyRFdp(double Vnorm, double f, double phase, const RMtrx& R1, bool full_acc)
            : Vn(Vnorm),k(twoPi*f/SpeedOfLight),phi0(phase),R(R1),fullacc(full_acc)
    {
        cosPhi0=cos(phi0);
        d0=1+Vn*cosPhi0;
    };

    void operator()(PSvector& p) const {
        R.Apply(p);
        if(fullacc)
            p.dp() = (p.dp()+Vn*(cos(phi0-k*p.ct())-cosPhi0))/d0;
        else
            p.dp() += Vn*cos(phi0-k*p.ct());
    }
};
};

//////////////////////////////////////////////////////////////////////////////
// Integrator classes
//////////////////////////////////////////////////////////////////////////////

namespace ParticleTracking{

namespace THIN_LENS {

// Integrator set definition
DEF_INTG_SET(ParticleComponentTracker,StdISet)
ADD_INTG(DriftCI)
ADD_INTG(SectorBendCI)
ADD_INTG(RectMultipoleCI)
ADD_INTG(LCAVIntegrator)
ADD_INTG(TransRFIntegrator)
ADD_INTG(MonitorCI)
ADD_INTG(SolenoidCI)
ADD_INTG(MarkerCI)
ADD_INTG(ParticleMapCI);
END_INTG_SET;

void DriftCI::TrackStep (double ds)
{
    CHK_ZERO(ds);
    ApplyDrift(currentBunch->GetParticles(),ds);
    return;
}

// Class TWRFStructureCI

void TWRFStructureCI::TrackStep (double ds)
{
    CHK_ZERO(ds);

    // Note that for particle tracking we use a higher order
    // approximation for the momentum error, and a linear matrix
    // for the transverse planes.

    const RFAcceleratingField& field = currentComponent->GetField();

    double g   = field.GetAmplitude();
    double f   = field.GetFrequency();
    double phi = field.GetPhase();

    if(g==0) {
        // cavity is off!
        ApplyDrift(currentBunch->GetParticles(),ds);
        return;
    }

    assert(f!=0);

    RMtrx Rm(2);
    double E0  = currentBunch->GetReferenceMomentum();

    TransportMatrix::TWRFCavity(ds,g,f,phi,E0,true,Rm.R);

    for_each(currentBunch->begin(),currentBunch->end(),ApplyRFdp(g*ds/E0,f,phi,Rm,true));

    if(true)
        currentBunch->IncrReferenceMomentum(g*ds*cos(phi));

    return;
}

// Class SectorBendCI

void SectorBendCI::TrackStep (double ds)
{
    CHK_ZERO(ds);

    double h = (*currentComponent).GetGeometry().GetCurvature();
    const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();

    double tilt = (*currentComponent).GetGeometry().GetTilt();

    if(tilt!=0) {
        RMtrx Rr;
        TransportMatrix::Srot(tilt,Rr.R);
        Rr.Apply(currentBunch->GetParticles());
    }

    //		if(GetIntegratedLength()==0 && pfi.entrance!=0)
    //			ApplyPoleFaceRotation(h,*pfi.entrance);

    MultipoleField& field = currentComponent->GetField();
    const double P0 = currentBunch->GetReferenceMomentum();
    const double q = currentBunch->GetChargeSign();
    const double Pref = currentComponent->GetMatchedMomentum(q);
    const double brho = P0/eV/SpeedOfLight;
    int np = field.HighestMultipole();

    assert(Pref>0);

    const Complex b0 = field.GetCoefficient(0);
    const Complex K1 = (np>0)? q*field.GetKn(1,brho) : Complex(0);

    // we need to split the magnet for a kick if the
    // following is true
    bool splitMagnet = b0.imag()!=0 || K1.imag()!=0 || np>1;
    double len = splitMagnet ? ds/2.0 : ds;

    RMtrx M(3,Pref);
    TransportMatrix::SectorBend(len,h,K1.real(),M.R);
    M.Apply(currentBunch->GetParticles(),P0);

    // Now if we have split the magnet, we need to
    // apply the kick approximation, and then
    // re-apply the map M
    if(splitMagnet) {

        // First we set the real parts of the
        // dipole and quad fields to zero, since these
        // components have been modeled in the matrix
        Complex b1=field.GetCoefficient(1);
        field.SetCoefficient(0,Complex(0,b0.imag()));
        field.SetCoefficient(1,Complex(0,b1.imag()));

        // Apply the integrated kick, and then track
        // through the linear second half
        for_each(currentBunch->begin(),currentBunch->end(),MultipoleKick(field,ds,P0,q));
        M.Apply(currentBunch->GetParticles(),P0);

        // Remember to set the components back
        field.SetCoefficient(0,b0);
        field.SetCoefficient(1,b1);
    }

    //		double Sr=IncrStep(ds);
    //		if(Sr==0 && pfi.exit!=0 )
    //			ApplyPoleFaceRotation(h,*pfi.exit);

    if(tilt!=0) {
        RMtrx Rr;
        TransportMatrix::Srot(-tilt,Rr.R);
        Rr.Apply(currentBunch->GetParticles());
    }

    return;
}

void SectorBendCI::TrackEntrance()
{
    double h = (*currentComponent).GetGeometry().GetCurvature();
    const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
    if(pfi.entrance!=0)
        ApplyPoleFaceRotation(h,*pfi.entrance);
}

void SectorBendCI::TrackExit()
{
    double h = (*currentComponent).GetGeometry().GetCurvature();
    const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
    if(pfi.exit!=0)
        ApplyPoleFaceRotation(h,*pfi.exit);
}

void SectorBendCI::ApplyPoleFaceRotation (double h, const SectorBend::PoleFace& pf)
{
    RMtrx M(2);
    TransportMatrix::PoleFaceRot(h,pf.rot,pf.fint,pf.hgap,M.R);
    M.Apply(currentBunch->GetParticles());
}

// Class RectMultipoleCI

void RectMultipoleCI::TrackStep (double ds)
{

    // Here we use a matrix to represent the quadrupole term, and a
    // single kick at the centre of the element for the other multipoles,
    // including any dipole term.
    using namespace TLAS;

    CHK_ZERO(ds);

    double P0 = currentBunch->GetReferenceMomentum();
    double q = currentBunch->GetChargeSign();
    double brho = P0/eV/SpeedOfLight;

    MultipoleField& field = currentComponent->GetField();

    const Complex cK1 = q*field.GetKn(1,brho);
    bool splitMagnet = field.GetCoefficient(0)!=0.0 || field.HighestMultipole()>1;
    double len = splitMagnet ? ds/2 : ds;

    RdpMtrx M(2);

    if(cK1!=0.0) {
        double K1 = cK1.imag()==0 ? cK1.real() : abs(cK1);
        TransportMatrix::QuadrupoleR(len,K1,M.R);
        TransportMatrix::QuadrupoleT(len,K1,M.T);

        if(cK1.imag()!=0) {
            // Need to rotate the map
            double a = arg(cK1)/2;
            RealMatrix Rr(4,4);
            TransportMatrix::Srot(a,Rr);
            M.R = Rr*M.R*Transpose(Rr);
            M.T = Rr*M.T*Transpose(Rr);
        }

        M.Apply(currentBunch->GetParticles());
    }
    else
        ApplyDrift(currentBunch->GetParticles(),len);

    if(splitMagnet) {
        Complex b1 = field.GetCoefficient(1);
        field.SetCoefficient(1,Complex(0));
        for_each(currentBunch->begin(),currentBunch->end(),MultipoleKick(field,ds,P0,q));
        if(b1!=0.0)
            M.Apply(currentBunch->GetParticles());
        else
            ApplyDrift(currentBunch->GetParticles(),len);
        field.SetCoefficient(1,b1);
    }
}

// Class SWRFStructureCI

void SWRFStructureCI::TrackStep (double ds)
{
    // Here we have a problem since the Chamber's matrix that we use
    // for the transport is only valid for n*(lambda/2), where n is
    // an integer. For now, we simple force ds to be a fixed number
    // of wavelengths

    CHK_ZERO(ds);

    const RFAcceleratingField& field = currentComponent->GetField();
    double g = field.GetAmplitude();
    double f = field.GetFrequency();
    double phi = field.GetPhase();
    double E0 = currentBunch->GetReferenceMomentum();

    double lambdaOver2 = SpeedOfLight/f/2;
    int ncells = static_cast<int>(ds/lambdaOver2);
    assert(ncells*lambdaOver2 == ds);

    RMtrx Rm(2);
    TransportMatrix::SWRFCavity(ncells,g,f,phi,E0,Rm.R);

    for_each(currentBunch->begin(),currentBunch->end(),ApplyRFdp(g*ds/E0,f,phi,Rm,true));

    if(true)
        currentBunch->IncrReferenceMomentum(g*ds*cos(phi));

    return;
}

// Class ExactRectMultipoleCI
/*****************************************************************
void ExactRectMultipoleCI::TrackStep (double ds)
{
	
	// Here we use a matrix to represent the quadrupole term, and a
	// single kick at the centre of the element for the other multipoles,
	// including any dipole term.
	using namespace TLAS;
	
	CHK_ZERO(ds);
	
	double P0 = GetBunch().GetReferenceMomentum();
	double q = GetBunch().GetChargeSign();
	double brho = P0/eV/SpeedOfLight;
	
	MultipoleField& field = currentComponent->GetField();
	
	const Complex cK1 = q*field.GetKn(1,brho);
	
	bool splitMagnet = field.GetCoefficient(0)!=0.0 || field.HighestMultipole()>1;
	double len = splitMagnet ? ds/2 : ds;
	
	if(cK1!=0.0) {
		double K1 = cK1.imag()==0 ? cK1.real() : abs(cK1);	
		// Now we iterate over all the particles in the bunch,
		// and calculate the exact map for each
		for(PSvectorArray::iterator p = GetBunch().begin(); p!=GetBunch().end(); p++) {
			RMtrx M(2);
			TransportMatrix::QuadrupoleR(len,K1/(1+p->dp()),M.R);			
			if(cK1.imag()!=0) {
				// Need to rotate the map
				double a = arg(cK1);
				RealMatrix Rr(4,4);
				TransportMatrix::Srot(a,Rr);
				M.R = Rr*M.R*Transpose(Rr);
			}
			M.Apply(*p);
		}
	}
	else 
		ApplyDrift(GetBunch().GetParticles(),len);
	
	
	if(splitMagnet) {
		
		Complex b1 = field.GetCoefficient(1);
		field.SetCoefficient(1,Complex(0));
		for_each(GetBunch().begin(),GetBunch().end(),MultipoleKick(field,ds,P0,q));
		field.SetCoefficient(1,b1);
		
		// Now repeat linear tracking through second half
		if(cK1!=0.0) {
			double K1 = cK1.imag()==0 ? cK1.real() : abs(cK1);	
			for(PSvectorArray::iterator p = GetBunch().begin(); p!=GetBunch().end(); p++) {
				RMtrx M(2);
				TransportMatrix::QuadrupoleR(len,K1/(1+p->dp()),M.R);			
				if(cK1.imag()!=0) {
					// Need to rotate the map
					double a = arg(cK1);
					RealMatrix Rr(4,4);
					TransportMatrix::Srot(a,Rr);
					M.R = Rr*M.R*Transpose(Rr);
				}
				M.Apply(*p);
			}
		}
		else 
			ApplyDrift(GetBunch().GetParticles(),len);	
	}
					  
	return;
}
*****************************************************************/

}; // end namespace THIN_LENS

// Class MonitorCI

void MonitorCI::TrackStep (double ds)
{

    double len = currentComponent->GetLength();
    if(len==0) {
        assert(ds==0);
        currentComponent->MakeMeasurement(*currentBunch);
    }
    else {
        double mpt=currentComponent->GetMeasurementPt()+len/2;
        if(GetIntegratedLength()+ds>mpt) {
            double s1=mpt-GetIntegratedLength();
            ApplyDrift(currentBunch->GetParticles(),s1);
            //				IncrStep(s1); // need to increment bunch timing
            currentComponent->MakeMeasurement(*currentBunch);
            ds-=s1;
        }
        ApplyDrift(currentBunch->GetParticles(),ds);
    }
    return;
}

// Class SolenoidCI

void SolenoidCI::TrackStep (double ds)
{
	const bool linear_map_only = false;

    CHK_ZERO(ds);
    double P0 = currentBunch->GetReferenceMomentum();
    double q = currentBunch->GetChargeSign();
    double brho = P0/eV/SpeedOfLight;

    double Bz = currentComponent->GetBz();

    if(fequal(Bz,0))
        ApplyDrift(currentBunch->GetParticles(),ds);
    else {
		if(linear_map_only) {
			RMtrx M(2);
			TransportMatrix::Solenoid(ds,q*Bz/brho,0,true,true,M.R);
			M.Apply(currentBunch->GetParticles());
		}
		else {
			// We use the exact momentum map for each particle energy.
			for(ParticleBunch::iterator p=currentBunch->begin(); p!=currentBunch->end(); p++) {
			RMtrx M(2);
			TransportMatrix::Solenoid(ds,q*Bz/brho/(1+p->dp()),0,true,true,M.R);
			M.Apply(*p);
			}
		}
    }
    return;
}

void MarkerCI::TrackStep(double)
{
    return;
}

}; // end namespace ParticleTracking


