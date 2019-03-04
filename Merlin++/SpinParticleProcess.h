/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: (c) 2004 Daniel A. Bates (LBNL) -- All Rights Reserved --
 */

#ifndef SpinParticleProcess_h
#define SpinParticleProcess_h 1

#include "ParticleBunch.h"
#include "ParticleBunchProcess.h"
#include "Components.h"

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

class SpinParticleBunch: public ParticleBunch
{
public:
	SpinParticleBunch(double P0, double Qm = 1);
	virtual ParticleBunch::iterator erase(ParticleBunch::iterator p);
	virtual size_t AddParticle(const Particle& p);
	size_t AddParticle(const Particle& p, const SpinVector& spin);
	virtual void push_back(const Particle& p);
	virtual void SortByCT();
	SpinVectorArray::iterator beginSpinArray();
	SpinVectorArray::iterator endSpinArray();
	virtual void Output(std::ostream& os) const;
	SpinVector GetAverageSpin() const;
	virtual bool ApplyTransformation(const Transform3D& t);

private:
	SpinVectorArray spinArray;
};

class SpinParticleProcess: public ParticleBunchProcess
{
public:
	SpinParticleProcess(int prio, int nstep = 1);

	/**
	 *	Sets the current accelerator component.
	 */
	virtual void SetCurrentComponent(AcceleratorComponent& component);

	/**
	 *	Preform the process for the specified step ds.
	 */
	virtual void DoProcess(double ds);

	/**
	 *	Returns the current maximum step length for this process.
	 *	@return Process current maximum step length
	 */
	virtual double GetMaxAllowedStepSize() const;

	/**
	 *	Sets the minimum number of equal steps to take through
	 *	the component.
	 */
	void SetNumComponentSteps(int n);

	/**
	 *  Set a momentum for calculating the spin precession
	 */
	void SetSpinMomentum(double p_spin);

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
