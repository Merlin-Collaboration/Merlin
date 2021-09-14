/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef EquilibriumDistribution_h
#define EquilibriumDistribution_h 1

// A.Wolski 19/11/01

#include "AcceleratorModel.h"
#include "SectorBend.h"

/**
 * Routines to calculate the equilibrium beam distribution
 * in an electron storage ring, using Chao's Method
 * J.Appl.Phys. 50(2), 1979
 *
 * A.Wolski, 28 June 2004.
 */

class IntegrateEigenvector
{
public:
	virtual ~IntegrateEigenvector()
	{
	}
	double Integral(ComplexVector& Ek, SectorBend* sb, double p0);

//protected:
	ComplexVector ek;
	Complex k;
	double h, tanE1;

//	void interpolate(double x[], double y[], int n, double X, double& Y, double& dy);
	double trapezium(double a, double b, int n);
	double romberg(double a, double b);
	virtual double function(double s) = 0;
};

class IntegrateZeroGradient: public IntegrateEigenvector
{
public:
	virtual ~IntegrateZeroGradient()
	{
	}
	IntegrateZeroGradient();
protected:
	virtual double function(double s);
};

class IntegrateWithGradient: public IntegrateEigenvector
{
public:
	virtual ~IntegrateWithGradient()
	{
	}
	IntegrateWithGradient();
protected:
	virtual double function(double s);
};

class EquilibriumDistribution
{
public:
	EquilibriumDistribution(AcceleratorModel* aModel, double refMomentum);
	double Emittance(int n);
	double DampingTime(int n);
	double DampingConstant(int n);
	double Tune(int n);
	double SynchronousTime();
	double BeamMoment(int i, int j, int ncpt = 0);

	void CalculateDampingConstants();
	void CalculateEmittance();

private:
	AcceleratorModel* theModel;
	double p0;
	double dampingConstant[3];
	double tune[3];
	double emittance[3];
	double synchronousTime;

};

#endif
