// (c) 2004 Daniel A. Bates (LBNL) -- All Rights Reserved --
// Part of the Spin Particle Process Package 1.0
// Modified 02/27/2004
// Further Modified A.Wolski 04/27/2004 (general tidying up)

#include "BeamDynamics/ParticleTracking/SpinParticleProcess.h"
#include "AcceleratorModel/StdComponent/SectorBend.h"
#include "AcceleratorModel/StdComponent/Solenoid.h"
#include "EuclideanGeometry/Space3D.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/utils.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

SpinVector::SpinVector()
{
    spin[0] = 0;
    spin[1] = 0;
    spin[2] = 1;
}

SpinVector::SpinVector(double _x, double _y, double _z)
{
    spin[0] = _x;
    spin[1] = _y;
    spin[2] = _z;
}

double SpinVector::x() const
{
    return spin[0];
}

double SpinVector::y() const
{
    return spin[1];
}

double SpinVector::z() const
{
    return spin[2];
}

double& SpinVector::x()
{
    return spin[0];
}

double& SpinVector::y()
{
    return spin[1];
}

double& SpinVector::z()
{
    return spin[2];
}

SpinParticleBunch::SpinParticleBunch(double P0, double Qm)
        : ParticleBunch(P0,Qm)
{}

size_t SpinParticleBunch::AddParticle(const Particle& p)
{
    SpinVector* spin = new SpinVector(0,0,1);
    spinArray.push_back(*spin);
    return ParticleBunch::AddParticle(p);
}

size_t SpinParticleBunch::AddParticle(const Particle& p, const SpinVector& spin)
{
    spinArray.push_back(spin);
    return ParticleBunch::AddParticle(p);
}

void SpinParticleBunch::push_back(const Particle& p)
{
    AddParticle(p);
}

void SpinParticleBunch::SortByCT()
{
    // Not yet fully implemented!
    // At present, only the list of phase space vectors is sorted;
    // the correlation with the list of spin vectors is lost.

    ParticleBunch::SortByCT();
}

void SpinParticleBunch::Output (std::ostream& os) const
{
    int oldp=os.precision(10);
    ios_base::fmtflags oflg = os.setf(ios::scientific,ios::floatfield);
	vector<SpinVector>::const_iterator sp = spinArray.begin();
    for(PSvectorArray::const_iterator p = begin(); p!=end(); p++,sp++) {
        os<<std::setw(24)<<GetReferenceTime();
        os<<std::setw(24)<<GetReferenceMomentum();
        for(size_t k=0; k<6; k++)
            os<<std::setw(20)<<(*p)[k];
		os<<std::setw(20)<<sp->x();
		os<<std::setw(20)<<sp->y();
		os<<std::setw(20)<<sp->z();
        os<<endl;
    }
    os.precision(oldp);
    os.flags(oflg);
}

SpinVector SpinParticleBunch::GetAverageSpin() const
{
	SpinVector pa(0,0,0);

    for(vector<SpinVector>::const_iterator sp = spinArray.begin(); sp!=spinArray.end(); sp++) {
		pa.x() += sp->x();
		pa.y() += sp->y();
		pa.z() += sp->z();
    }
	pa.x() /= spinArray.size();
	pa.y() /= spinArray.size();
	pa.z() /= spinArray.size();

	return pa;
}

ParticleBunch::iterator SpinParticleBunch::erase(ParticleBunch::iterator p)
{
    // find the 'offset' of p from the start of the particle bunch
    size_t n = distance(pArray.begin(),p);

    // remove the n-th spin vector
    SpinVectorArray::iterator it = spinArray.begin();
    advance(it,n);
    spinArray.erase(it);

    // called the base function to remove the PSvector
    return ParticleBunch::erase(p);
}

bool SpinParticleBunch::ApplyTransformation (const Transform3D& t)
{
	// First apply the transformation to the particle coordinates
	ParticleBunch::ApplyTransformation(t);

	// If the transformation represents a rotation, then we must also
	// rotate the spin vectors. (Note that we must now iterate twice
	// through the particles which is a little inefficient.)
	if(!t.R().isIdentity()) {
		const Rotation3D& R(t.R());
		for(vector<SpinVector>::iterator sp = spinArray.begin(); sp!=spinArray.end(); sp++) {
			Vector3D S(sp->x(),sp->y(),sp->z());
			S=R(S);
			sp->x()=S.x;
			sp->y()=S.y;
			sp->z()=S.z;
		}
	}
	return true;
}


vector<SpinVector>::iterator SpinParticleBunch::beginSpinArray()
{
    return spinArray.begin();
}

vector<SpinVector>::iterator SpinParticleBunch::endSpinArray()
{
    return spinArray.end();
}

SpinParticleProcess::SpinParticleProcess(int prio, int nstep)
        : ParticleBunchProcess("SPIN TRACKING PROCESS",prio),ns(nstep),pspin(0) {}

void SpinParticleProcess::SetCurrentComponent(AcceleratorComponent& component){
    nk1 = 0;
    intS = 0;

    //Get current field of current component
    currentField = component.GetEMField();

    //Get length of current component/ns
    clength = component.GetLength();
    dL = clength/ns;

    //Determine if the component is a SectorBend
    sbend = dynamic_cast<SectorBend*>(&component);

    //Determine if the component is a Solenoid
    solnd = dynamic_cast<Solenoid*>(&component);

    //Determine if the present bunch has any spin information
    SpinParticleBunch* spinbunch = dynamic_cast<SpinParticleBunch*>(currentBunch);

    if (currentField && spinbunch) {
        active = true;
    } else {
        active = false;
    }
}

struct RotateSpinVector
{
private:
    bool isBend;

public:
    RotateSpinVector(bool bend)
            : isBend(bend)
    {};

    void RotateSpin(const Vector3D& bnorm, double ds, SpinVector& spin, double gamma)
    {
        double& spinx = spin.x();
        double& spiny = spin.y();
        double& spinz = spin.z();

        Vector3D w;
        if (isBend) {	// Arc geometry
            w.x = -(1+ElectronGe*gamma)*bnorm.x;
            w.y = -(  ElectronGe*gamma)*bnorm.y;
            w.z = -(1+ElectronGe      )*bnorm.z;
        } else {		// Rectangular geometry
            w.x = -(1+ElectronGe*gamma)*bnorm.x;
            w.y = -(1+ElectronGe*gamma)*bnorm.y;
            w.z = -(1+ElectronGe      )*bnorm.z;
        }

        // Calculate Cross Products
        double pX_crossproduct = w.y*spinz - w.z*spiny;
        double pY_crossproduct = w.z*spinx - w.x*spinz;
        double pZ_crossproduct = w.x*spiny - w.y*spinx;

        // Calculate Scalar Product
        double scalarproduct = w.z*spinz + w.y*spiny + w.x*spinx;

        // Calculate and apply spin rotation
        double w2    = w.x*w.x + w.y*w.y + w.z*w.z;
        double omega = sqrt(w2);

        double dt = ds / SpeedOfLight;
        double cosOmegaT = cos(dt*omega);
        double sinOmegaT = sin(dt*omega);

        double pX = spinx*cosOmegaT + w.x*scalarproduct*(1-cosOmegaT)/w2 + pX_crossproduct*sinOmegaT/omega;
        double pY = spiny*cosOmegaT + w.y*scalarproduct*(1-cosOmegaT)/w2 + pY_crossproduct*sinOmegaT/omega;
        double pZ = spinz*cosOmegaT + w.z*scalarproduct*(1-cosOmegaT)/w2 + pZ_crossproduct*sinOmegaT/omega;

        spinx = pX;
        spiny = pY;
        spinz = pZ;
    };

};

void SpinParticleProcess::SetSpinMomentum(double p_spin)
{
    pspin = p_spin;
};

void SpinParticleProcess::DoProcess(double ds)
{
    double P0   = currentBunch->GetReferenceMomentum();
    double brho = P0/eV/SpeedOfLight;
    bool isBend = sbend ? true : false;
	bool isSolenoid = solnd ? true : false;

    if(pspin != 0)
        P0 = pspin;

    RotateSpinVector rot(isBend);

    Vector3D b;

    SpinParticleBunch* spinbunch = dynamic_cast<SpinParticleBunch*>(currentBunch);
    SpinVectorArray::iterator spin = spinbunch->beginSpinArray();
    for(PSvectorArray::iterator p = spinbunch->begin(); p!=spinbunch->end(); p++,spin++) {
        double norm  = SpeedOfLight/brho/(1.0+p->dp());
        double gamma = P0*(1.0+p->dp())/(ElectronMassMeV*MeV);

		if (intS==0) {

			// Apply spin rotation from dipole entrance fringe field
			if (isBend) {
				SectorBend::PoleFace* pf = sbend->GetPoleFaceInfo().entrance;
				double theta = pf ? pf->rot : 0;
				double intbz = 0.5*sbend->GetB0()*p->y()*norm;
				b.x = sin(theta)*intbz;
				b.y = 0;
				b.z = cos(theta)*intbz;
				if ( b.x!=0 || b.z!=0 )
					rot.RotateSpin(b, 1.0, *spin, gamma);
			}

			// Apply spin rotation from solenoid entrance fringe field
			// We use a hard-edged model for the fringe field;
			// a positive value for the solenoid field means the field
			// is pointing in the direction of the beam.
			if (isSolenoid) {
				double bz = solnd->GetBz();
				b.x = -bz*p->x();
				b.y = -bz*p->y();
				b.z = 0.0;
				if ( b.x!=0 || b.y!=0 )
					rot.RotateSpin(b, 1.0, *spin, gamma);
			}
		}

        // Apply spin rotation from body of magnet
        b = currentField->GetBFieldAt(Point3D(p->x(),p->y(),0));

        if ( b.x!=0 || b.y!=0 || b.z!=0 ) {
            b.x *= norm;
            b.y *= norm;
            b.z *= norm;
            rot.RotateSpin(b, ds, *spin, gamma);
        }

		if (fequal(intS+ds,clength)) {

			// Apply spin rotation from dipole exit fringe field
			if (isBend) {
				SectorBend::PoleFace* pf = sbend->GetPoleFaceInfo().exit;
				double theta = pf ? pf->rot : 0;
				double intbz = 0.5*sbend->GetB0()*p->y()*norm;
				b.x = sin(theta)*intbz;
				b.y = 0;
				b.z =-cos(theta)*intbz;
				if ( b.x!=0 || b.z!=0 )
					rot.RotateSpin(b, 1.0, *spin, gamma);
			}

			// Apply spin rotation from solenoid exit fringe field
			// We use a hard-edged model for the fringe field
			if (isSolenoid) {
				double bz = solnd->GetBz();
				b.x = bz*p->x();
				b.y = bz*p->y();
				b.z = 0.0;
				if ( b.x!=0 || b.y!=0 )
					rot.RotateSpin(b, 1.0, *spin, gamma);
			}
		}

    };

    intS += ds;

}

double SpinParticleProcess::GetMaxAllowedStepSize() const
{
    return (nk1+1)*dL-intS;
}

void SpinParticleProcess::SetNumComponentSteps(int n)
{
    ns = n;
}
