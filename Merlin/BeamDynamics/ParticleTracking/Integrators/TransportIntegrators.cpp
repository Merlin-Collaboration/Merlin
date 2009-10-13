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
// $Revision: 1.11 $
// 
/////////////////////////////////////////////////////////////////////////

#include "BeamDynamics/ParticleTracking/Integrators/StdIntegrators.h"
#include "BeamDynamics/ParticleTracking/Integrators/LCAVintegrator.h"
#include "BasicTransport/BasicTransportMaps.h"
#include "BasicTransport/TransportMatrix.h"
#include "BasicTransport/MatrixMaps.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "BeamDynamics/ParticleTracking/Integrators/TransRFIntegrator.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

namespace ParticleTracking {
namespace TRANSPORT {

// Integrator set definition
DEF_INTG_SET(ParticleComponentTracker,StdISet)
ADD_INTG(DriftCI)
ADD_INTG(SectorBendCI)
ADD_INTG(RectMultipoleCI)
ADD_INTG(LCAVIntegrator)
ADD_INTG(TransRFIntegrator)
ADD_INTG(SolenoidCI)
ADD_INTG(THIN_LENS::SWRFStructureCI)
ADD_INTG(MarkerCI)
ADD_INTG(MonitorCI)
ADD_INTG(ParticleMapCI);
END_INTG_SET;

}; // end namespace TRANSPORT
}; // end namespace ParticleTracking


template<> MAKE_DEF_INTG_SET(ParticleTracking::ParticleComponentTracker,ParticleTracking::TRANSPORT::StdISet);

#define CHK_ZERO(s) if(s==0) return;

namespace ParticleTracking {

namespace {

// tolerance for bend scaling
#define REL_ENGY_TOL 1.0-06

struct MultipoleKick {
    const MultipoleField& field;
    Complex scale;

    MultipoleKick(const MultipoleField& f, double len, double P0, double q, double phi =0)
            : field(f)
    {
        scale = q*len*eV*SpeedOfLight/P0*Complex(cos(phi),sin(phi));
    }

    void operator()(PSvector& v) {
        double x=v.x();
        double y=v.y();
        double dp=v.dp();
        Complex F = scale*field.GetField2D(x,y)/(1+dp);
        v.xp() += -F.real();
        v.yp() +=  F.imag();
    };
};

// Apply a map without dp/p scaling
struct ApplyMap {
    RTMap* m;
    ApplyMap(RTMap* amap) : m(amap) {}
    void operator()(PSvector& p) {
        m->Apply(p);
    }
};

// Apply map with a dp/p scaling
struct ApplyMap1 {
    RTMap* m;
    double Eratio;

    ApplyMap1(RTMap* amap, double Er) : m(amap), Eratio(Er) {}
    void operator()(PSvector& p) {
        double dp = p.dp();
        p.dp() = Eratio*(1+dp)-1;
        m->Apply(p);
        p.dp() = dp;
    }
};

struct ApplyDrift {
    const double s;
    ApplyDrift(double len) : s(len) {}
    void operator()(PSvector& p) {
        const double xp = p.xp();
        const double yp = p.yp();
        p.x()+=s*xp;
        p.y()+=s*yp;
        p.ct()-=s*(xp*xp+yp*yp)/2.0;
    }
};

// Functor ApplyRFdp (used for full acceleration)
struct ApplyRFMap {

    double Vn,k,phi0,cosPhi0,d0;
    RTMap* m;
    bool fullacc;

    ApplyRFMap(double Vnorm, double f, double phase, RTMap* m1, bool full_acc)
            : Vn(Vnorm),k(twoPi*f/SpeedOfLight),phi0(phase),m(m1),fullacc(full_acc)
    {
        cosPhi0=cos(phi0);
        d0=1+Vn*cosPhi0;
    };

    void operator()(PSvector& p) const {
        m->Apply(p);
        if(fullacc)
            p.dp() = (p.dp()+Vn*(cos(phi0-k*p.ct())-cosPhi0))/d0;
        else
            p.dp() += Vn*cos(phi0-k*p.ct());
    }
};


inline void ApplyMapToBunch(ParticleBunch& bunch, RTMap* amap)
{
    for_each(bunch.begin(),bunch.end(),ApplyMap(amap));
}

inline void ApplyMapToBunch(ParticleBunch& bunch, RTMap* amap, double Er)
{
    for_each(bunch.begin(),bunch.end(),ApplyMap1(amap,Er));
}

inline void ApplyDriftToBunch(ParticleBunch& bunch, double len)
{
    for_each(bunch.begin(),bunch.end(),ApplyDrift(len));
}

void RotateBunchAboutZ(ParticleBunch& bunch, double phi)
{
    RMtrx M(2);
    TransportMatrix::Srot(phi,M.R);
    M.Apply(bunch.GetParticles());
}

inline bool operator==(const Complex& z, double x)
{
    return z.imag()==0 && z.real()==x;
}

inline bool operator!=(const Complex& z, double x)
{
    return !(z==x);
}
};

namespace TRANSPORT {

void DriftCI::TrackStep (double ds)
{
    CHK_ZERO(ds);
    RTMap* m = DriftTM(ds);
    ApplyMapToBunch(*currentBunch,m);
    delete m;
    return;
}


void SectorBendCI::TrackStep(double ds)
{
    CHK_ZERO(ds);

    double h = (*currentComponent).GetGeometry().GetCurvature();
    MultipoleField& field = (*currentComponent).GetField();
    const double P0 = (*currentBunch).GetReferenceMomentum();
    const double q = (*currentBunch).GetChargeSign();
    const double Pref = (*currentComponent).GetMatchedMomentum(q);
    const double brho = P0/eV/SpeedOfLight;
    int np = field.HighestMultipole();

    assert(Pref>0);

    const Complex b0 = field.GetCoefficient(0);
    const Complex K1 = (np>0)? q*field.GetKn(1,brho) : Complex(0);

    // we need to split the magnet for a kick if the
    // following is true
    bool splitMagnet = b0.imag()!=0 || K1.imag()!=0 || np>1;
    double len = splitMagnet ? ds/2.0 : ds;

    // Construct the second-order map
    RTMap* M = (abs(K1)==0) ? SectorBendTM(len,h) : GenSectorBendTM(len,h,K1.real(),0);

    if(fequal(P0,Pref,REL_ENGY_TOL))
        ApplyMapToBunch(*currentBunch,M);
    else
        ApplyMapToBunch(*currentBunch,M,P0/Pref);

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
        for_each((*currentBunch).begin(),(*currentBunch).end(),MultipoleKick(field,ds,P0,q));

        if(fequal(P0,Pref,REL_ENGY_TOL))
            ApplyMapToBunch(*currentBunch,M);
        else
            ApplyMapToBunch(*currentBunch,M,P0/Pref);

        // Remember to set the components back
        field.SetCoefficient(0,b0);
        field.SetCoefficient(1,b1);
    }

    // must delete the map
    delete M;
    return;

}

void SectorBendCI::TrackEntrance()
{
    const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
    double tilt = (*currentComponent).GetGeometry().GetTilt();
    if(tilt!=0) {
        //	cout<<"rotating by "<<tilt<<endl;
        RotateBunchAboutZ(*currentBunch,-tilt);
    }
    ApplyPoleFaceRotation(pfi.entrance);
}

void SectorBendCI::TrackExit()
{
    const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
    double tilt = (*currentComponent).GetGeometry().GetTilt();
    ApplyPoleFaceRotation(pfi.exit);
    if(tilt!=0){
        //	cout<<"rotating by "<<-tilt<<endl;
        RotateBunchAboutZ(*currentBunch,tilt);
    }
}

void SectorBendCI::ApplyPoleFaceRotation (const SectorBend::PoleFace* pf)
{
#define _PFV(p,v) !(p) ? 0 : p->v;

    const double P0 = currentBunch->GetReferenceMomentum();
    const double brho = P0/eV/SpeedOfLight;
    const double h = (*currentComponent).GetGeometry().GetCurvature();
    const double k = (*currentComponent).GetField().GetKn(1,brho).real();

    double beta = _PFV(pf,rot);
    double c = 0; // currently not implemented
    double hg = _PFV(pf,hgap);
    double fint = _PFV(pf,fint);
    RTMap* M = PoleFaceTM(h,k,beta,c,fint,hg);
    ApplyMapToBunch(*currentBunch,M);
    delete M;
}


void RectMultipoleCI::TrackStep (double ds)
{
    // Here we use a matrix to represent the quadrupole term, and a
    // single kick at the centre of the element for the other multipoles,
    // including any dipole term.
    using namespace TLAS;



    double P0 = (*currentBunch).GetReferenceMomentum();
    double q = (*currentBunch).GetChargeSign();
    double brho = P0/eV/SpeedOfLight;
    MultipoleField& field = (*currentComponent).GetField();

    // we now support thin-lens kicks (this has been added to support
    // thin-lens corrector dipoles)

    if((*currentComponent).GetLength()==0 && ds==0 && !field.IsNullField()) {
        // treat field as integrated strength
        for_each((*currentBunch).begin(),(*currentBunch).end(),MultipoleKick(field,1.0,P0,q));
        return;
    }

    CHK_ZERO(ds);

    const Complex ch  = q*field.GetKn(0,brho);
    const Complex cK1 = q*field.GetKn(1,brho);
    const Complex cK2 = q*field.GetKn(2,brho);
    int np = field.HighestMultipole();

    // We split the magnet for a multipole kick if:
    //  - there is a dipole field
    //  - there is a quad field with higher-order multipoles
    //  - there is a sextupole    "     "     "        "

    bool splitMagnet = abs(ch)!=0 || (cK1!=Complex(0) && np>1) || np>2;
    double len = splitMagnet ? ds/2 : ds;


    if(cK1!=0.0) { // quad R+T matrix with thin-lens kicks for other multipoles
        double K1,phi;

        if(imag(cK1)==0) {
            K1 = real(cK1);
            phi=0;
        }
        else {
            K1 = abs(cK1);
            phi = arg(cK1)/2.0;
        }

        if(!fequal(phi,0))
            RotateBunchAboutZ(*currentBunch,-phi);

        RTMap* M = QuadrupoleTM(len,K1);
        ApplyMapToBunch(*currentBunch,M);
        if(splitMagnet) {
            Complex b1 = field.GetCoefficient(1);
            field.SetCoefficient(1,Complex(0));
            for_each((*currentBunch).begin(),(*currentBunch).end(),MultipoleKick(field,ds,P0,q,-phi));
            // Apply second half of map
            ApplyMapToBunch(*currentBunch,M);
            field.SetCoefficient(1,b1);
        }
        delete M;
        if(!fequal(phi,0))
            RotateBunchAboutZ(*currentBunch,phi);
    }
    else if(cK2!=0.0) { // sextupole R+T matrix with thin-lens kicks for other multipoles

        double K2,phi;
        if(imag(cK2)==0) {
            K2 = real(cK2);
            phi=0;
        }
        else {
            K2 = abs(cK2);
            phi = arg(cK2)/3.0;
        }

        if(!fequal(phi,0))
            RotateBunchAboutZ(*currentBunch,-phi);

        RTMap* M = SextupoleTM(len,K2);
        ApplyMapToBunch(*currentBunch,M);
        if(splitMagnet) {
            Complex b2 = field.GetCoefficient(2);
            field.SetCoefficient(2,Complex(0));
            for_each((*currentBunch).begin(),(*currentBunch).end(),MultipoleKick(field,ds,P0,q,-phi));
            // Apply second half of map
            ApplyMapToBunch(*currentBunch,M);
            field.SetCoefficient(2,b2);
        }
        delete M;
        if(!fequal(phi,0))
            RotateBunchAboutZ(*currentBunch,phi);
    }
    else { // drift with a kick in the middle
        ApplyDriftToBunch(*currentBunch,len);
        if(splitMagnet) {
            for_each((*currentBunch).begin(),(*currentBunch).end(),MultipoleKick(field,ds,P0,q));
            // Apply second half of map
            ApplyDriftToBunch(*currentBunch,len);
        }
    }

    return;
}

/*****
void SWRFStructureCI::TrackStep(double ds)
{
	
	// Here we have a problem since the Chamber's matrix that we use
	// for the transport is only valid for n*(lambda/2), where n is 
	// an integer. For now, we simple force ds to be a fixed number
	// of wavelengths
	
	CHK_ZERO(ds);
	
	const SWRFfield& field = (*currentComponent).GetField();
	double g = field.GetAmplitude();
	double f = field.GetFrequency();
	double phi = field.GetPhase();
	double E0 = (*currentBunch).GetReferenceMomentum();
	
	double lambdaOver2 = SpeedOfLight/f/2;
	int ncells = static_cast<int>(ds/lambdaOver2);
	assert(ncells*lambdaOver2 == ds);
	
	RealMatrix R(2,2);
	TransportMatrix::SWRFCavity(ncells,g,f,phi,E0,R);
	RTMap* M = new RTMap();
	
	(*M)(3,3)=(*M)(1,1)=R(0,0);
	(*M)(3,4)=(*M)(1,2)=R(0,1);
	(*M)(4,3)=(*M)(2,1)=R(1,0);
	(*M)(4,4)=(*M)(2,2)=R(1,1);
	(*M)(5,5)=(*M)(6,6)=1.0;
	
	for_each((*currentBunch).begin(),(*currentBunch).end(),ApplyRFMap(g*ds/E0,f,phi,M,true));
	
	if(true)
		(*currentBunch).IncrReferenceMomentum(g*ds*cos(phi));
	
	return;
}
**/

}; // end of namespace TRANSPORT
}; // end of namespace ParticleTracking

