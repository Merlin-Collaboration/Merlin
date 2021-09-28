/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <fstream>
#include "ParticleBunch.h"
#include "ParticleTracker.h"
#include "SynchRadParticleProcess.h"
#include "TransferMatrix.h"
#include "EquilibriumDistribution.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"
#include "MatrixPrinter.h"
#include "Interpolation.h"

// The following constants specify the integration precision
#define EPS 1.0e-8
#define MAXPOINTS 30
#define K 5

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;
using namespace ParticleTracking;

/*
 * void IntegrateEigenvector::interpolate(double* x, double* y, int n, double X, double& Y, double& dy)
   {
    int j = 1;
    double close, closest, ho, hp, w;
    double c[MAXPOINTS];
    double d[MAXPOINTS];

    closest = fabs(X - x[1]);

    for(int i = 1; i <= n; i++)
    {
        c[i] = d[i] = y[i];
        if((close = fabs(X - x[i])) < closest)
        {
            j = i;
            closest = close;
        }
    }

    Y = y[j--];

    for(int m = 1; m < n; m++)
    {
        for(int i = 1; i <= n; i++)
        {
            ho = x[i] - X;
            hp = x[i + m] - X;
            w = c[i + 1] - d[i];
            d[i] = hp * w/(ho - hp);
            c[i] = ho * w/(ho - hp);
        }
        Y += (dy = (2 * j < (n - m) ? c[j + 1] : d[j--]));
    }
   }
 */

double IntegrateEigenvector::trapezium(double a, double b, int n)
{
	double x, sum, delta;
	static double answer; // should not be modified between successive calls
	int npts;

	if(n == 1)
	{
		return answer = (b - a) * (function(a) + function(b)) / 2.0;
	}
	else
	{
		npts = 1 << (n - 2);
		delta = (b - a) / npts;
		x = a + 0.5 * delta;
		sum = 0;
		for(int j = 0; j < npts; j++, x += delta)
		{
			sum += function(x);
		}
		return answer = (answer + delta * sum) / 2.0;
	}
}

double IntegrateEigenvector::romberg(double a, double b)
{
	double value, error;
	double s[MAXPOINTS + 1];
	double h[MAXPOINTS + 2];
	int j;

	h[1] = 1.0;
	for(j = 1; j <= MAXPOINTS; j++)
	{
		s[j] = trapezium(a, b, j);
		if(j >= K)
		{
			Interpolation I(h + j - K, s + j - K, K, K);
			value = I.itsMethod->ValueAt(0.0, &error);
			if(fabs(error) <= EPS * fabs(value))
			{
				return value;
			}
		}
		h[j + 1] = 0.25 * h[j];
	}
	return 0;
}

double IntegrateEigenvector::Integral(ComplexVector& Ek, SectorBend* sb, double p0)
{
	ek = Ek;
	h = sb->GetB0() * eV * SpeedOfLight / p0;
	k = sqrt(Complex(sb->GetB1() * eV * SpeedOfLight / p0, 0));

	SectorBend::PoleFace* pf = (sb->GetPoleFaceInfo()).entrance;
	tanE1 = pf ? tan(pf->rot) : 0.0;

	return romberg(0, sb->GetLength());
}

IntegrateWithGradient::IntegrateWithGradient()
{
}

double IntegrateWithGradient::function(double s)
{
	Complex sbendTM52 = -(1.0 - cos(k * s)) * h / k / k;
	Complex sbendTM51 = -sin(k * s) * h / k + h * tanE1 * sbendTM52;
	Complex sbendTM56 = -pow(h / k, 2) * (s - sin(k * s) / k);

	Complex Ek5l = sbendTM51 * ek(0) + sbendTM52 * ek(1) + ek(4) + sbendTM56 * ek(5);
	return abs(Ek5l) * abs(Ek5l);
}

IntegrateZeroGradient::IntegrateZeroGradient()
{
}

double IntegrateZeroGradient::function(double s)
{
	Complex sbendTM52 = -s * s * h / 2.;
	Complex sbendTM51 = -s * h + h * tanE1 * sbendTM52;
	Complex sbendTM56 = sbendTM52 * s * h / 3.;

	Complex Ek5l = sbendTM51 * ek(0) + sbendTM52 * ek(1) + ek(4) + sbendTM56 * ek(5);
	return abs(Ek5l) * abs(Ek5l);
}

EquilibriumDistribution::EquilibriumDistribution(AcceleratorModel* aModel, double refMomentum) :
	theModel(aModel), p0(refMomentum)
{
	for(int n = 0; n < 3; n++)
	{
		dampingConstant[n] = 0;
		emittance[n] = 0;
		tune[n] = 0;
	}

	synchronousTime = 0;
}

double EquilibriumDistribution::DampingConstant(int n)
{
	return dampingConstant[n];
}

double EquilibriumDistribution::DampingTime(int n)
{
	double trev = (*theModel).GetGlobalFrame().GetGeometryLength() / SpeedOfLight;
	return trev / dampingConstant[n];
}

double EquilibriumDistribution::Tune(int n)
{
	return tune[n];
}

double EquilibriumDistribution::Emittance(int n)
{
	return emittance[n];
}

double EquilibriumDistribution::SynchronousTime()
{
	return synchronousTime;
}

void EquilibriumDistribution::CalculateDampingConstants()
{
	TransferMatrix tm(theModel, p0);
	tm.Radiation(true);
	tm.SetRadStepSize(0.1);

	RealMatrix M(6);
	PSvector p(0);
	tm.FindClosedOrbitTM(M, p);
	synchronousTime = p.ct() / SpeedOfLight;

	ComplexVector eigenvalues(3);
	ComplexMatrix eigenvectors(6, 3);
	EigenSystem(M, eigenvalues, eigenvectors);

	for(int m = 0; m < 3; m++)
	{
		tune[m] = arg(eigenvalues[m]) / twoPi;
		if(tune[m] < 0)
		{
			tune[m] = 1.0 + tune[m];
		}
		dampingConstant[m] = -log(abs(eigenvalues[m]));
	}
}

void EquilibriumDistribution::CalculateEmittance()
{
	const double dscale = 1.0e-7;

	int k; // used twice, but MS doesn't following ANSI scoping rules

	double SumE5[3];
	for(k = 0; k < 3; k++)
	{
		SumE5[k] = 0.0;
	}

	PSvector orbit(0);
	RealMatrix M(6);
	TransferMatrix tm(theModel, p0);
	tm.FindClosedOrbitTM(M, orbit);

	ofstream rmatrix("matrix.dat");
	rmatrix << "Before symplectic adjustment\n\n";
	MatrixForm(M, rmatrix, OPFormat().precision(12).fixed());

	// Need to fix the R53 and R54 terms, which are not calculated correctly
	// for a distorted closed orbit (tracking problem!)
	M(4, 2) = -(-M(2, 2) * M(3, 1) * M(4, 0) + M(2, 1) * M(3, 2) * M(4, 0)
		+ M(2, 2) * M(3, 0) * M(4, 1) - M(2, 0) * M(3, 2) * M(4, 1)
		+ M(2, 5) * M(3, 2) * M(4, 4) - M(2, 2) * M(3, 5) * M(4, 4)
		- M(2, 4) * M(3, 2) * M(4, 5) + M(2, 2) * M(3, 4) * M(4, 5))
		/ (M(2, 3) * M(3, 2) - M(2, 2) * M(3, 3));

	M(4, 3) = (-M(3, 3) * M(2, 1) * M(4, 0) + M(3, 3) * M(2, 0) * M(4, 1)
		- M(3, 3) * M(2, 5) * M(4, 4) + M(3, 3) * M(2, 4) * M(4, 5)
		+ M(2, 3) * M(3, 1) * M(4, 0) - M(2, 3) * M(3, 0) * M(4, 1)
		+ M(2, 3) * M(3, 5) * M(4, 4) - M(2, 3) * M(3, 4) * M(4, 5))
		/ (M(2, 3) * M(3, 2) - M(2, 2) * M(3, 3));

	//Just make sure the matrix is symplectic
	Symplectify(M);
	rmatrix << "\n\nAfter symplectic adjustment\n\n" << endl;
	MatrixForm(M, rmatrix, OPFormat().precision(12).fixed());

	Symplectify(M); //Just make sure the matrix is symplectic

	ComplexVector eigenvalues(3);
	ComplexMatrix eigenvectors(3, 6);
	ComplexVector newEigenvector(6);

	EigenSystem(M, eigenvalues, eigenvectors);

	ParticleBunch particle(p0, 1.0);
	for(int n = 0; n < 7; n++)
	{
		Particle p = orbit;
		if(n > 0)
		{
			p[n - 1] += dscale;
		}
		particle.push_back(p);
	}

	ParticleTracker tracker(theModel->GetBeamline(), &particle);
	tracker.InitStepper();
	bool loop = true;

	do
	{

		SectorBend* sb = dynamic_cast<SectorBend*>(&(tracker.GetCurrentComponent()));

		if(sb && sb->GetB0())
		{

			ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
			const Particle& p_ref = *ip++;

			for(k = 0; k < 3; k++)
			{

				newEigenvector = 0;

				for(int row = 0; row < 6; row++)
				{
					ip = tracker.GetTrackedBunch().begin();
					ip++;
					for(int col = 0; col < 6; col++, ip++)
					{
						newEigenvector(row) += ((*ip)[row] - p_ref[row]) * eigenvectors[k](col) / dscale;
					}
				}

				double intgrl = 0.0;

				if(sb->GetB1() == 0)
				{
					IntegrateZeroGradient integrator;
					intgrl = integrator.Integral(newEigenvector, sb, p0);
				}
				else
				{
					IntegrateWithGradient integrator;
					intgrl = integrator.Integral(newEigenvector, sb, p0);
				}

				SumE5[k] += intgrl / fabs(pow(p0 / eV / SpeedOfLight / sb->GetB0(), 3));
			}
		}

		loop = tracker.StepComponent();

	} while(loop);

	const double C_L = 55 * ElectronRadius * PlanckConstant / (48 * twoPi * sqrt(3.0) * ElectronMass);
	double gamma5 = pow(p0 * ElectronCharge / eV / ElectronMass / SpeedOfLight / SpeedOfLight, 5);

	for(k = 0; k < 3; k++)
	{
		emittance[k] = SumE5[k] * C_L * gamma5 / SpeedOfLight / dampingConstant[k];
	}
}

double EquilibriumDistribution::BeamMoment(int i, int j, int ncpt)
{
	ComplexVector eigenvalues(3);
	ComplexMatrix eigenvectors(6, 3);
	RealMatrix M(6);

	TransferMatrix tm(theModel, p0);
	tm.SetObservationPoint(ncpt);
	tm.FindTM(M);

	EigenSystem(M, eigenvalues, eigenvectors);

	double moment = 0;
	for(int k = 0; k < 3; k++)
		moment += 2 * emittance[k] * (eigenvectors(k, i).real() * eigenvectors(k, j).real()
			+ eigenvectors(k, i).imag() * eigenvectors(k, j).imag());

	return moment;
}
