// (c) 2004 Daniel A. Bates (LBNL) -- All Rights Reserved --
// Part of the Spin Particle Process Package 1.0
// Modified 02/27/2004
// Further Modified A.Wolski 04/27/2004 (general tidying up)

#ifndef SpinParticleProcess_h
#define SpinParticleProcess_h 1

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "AcceleratorModel/Components.h"

using namespace ParticleTracking;

class SpinVector
{
public:
    SpinVector();
    SpinVector(double _x, double _y, double _z);
    double x() const;
    double y() const;
    double z() const;
    double& x();
    double& y();
    double& z();

private:
    double spin[3];
};

typedef vector<SpinVector> SpinVectorArray;

class SpinParticleBunch : public ParticleBunch
{
public:
    SpinParticleBunch(double P0, double Qm=1);
    virtual ParticleBunch::iterator erase(ParticleBunch::iterator p);
    virtual size_t AddParticle(const Particle& p);
    size_t AddParticle(const Particle& p, const SpinVector& spin);
    virtual void push_back(const Particle& p);
    virtual void SortByCT();
    SpinVectorArray::iterator beginSpinArray();
    SpinVectorArray::iterator endSpinArray();
	virtual void Output (std::ostream& os) const;
	SpinVector GetAverageSpin() const;
    virtual bool ApplyTransformation (const Transform3D& t);

private:
    SpinVectorArray spinArray;
};


class SpinParticleProcess : public ParticleBunchProcess
{
public:
    SpinParticleProcess (int prio, int nstep = 1);
    //	Sets the current accelerator component.
    virtual void SetCurrentComponent (AcceleratorComponent& component);
    //	Preform the process for the specified step ds.
    virtual void DoProcess (double ds);
    //	Returns the current maximum step length for this process.
    virtual double GetMaxAllowedStepSize () const;
    //	Sets the minimum number of equal steps to take through
    //	the component.
    void SetNumComponentSteps (int n);
    //  Set a momentum for calculating the spin precession
    void SetSpinMomentum (double p_spin);

private:
    int ns;
    int nk1;
    double dL;
    double intS;
    const SectorBend* sbend;
	const Solenoid* solnd;
    const EMField* currentField;
    double clength;
    double pspin;
};

#endif
