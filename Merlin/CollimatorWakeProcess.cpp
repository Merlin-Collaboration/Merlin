/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <vector>

#include "merlin_config.h"

#include "Aperture.h"

#include "ParticleBunch.h"
#include "ParticleBunchUtilities.h"
#include "WakeFieldProcess.h"

#include "CollimatorWakePotentials.h"
#include "CollimatorWakeProcess.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

namespace
{
inline double powd(double x, int y)
{
	double product(1);
	for(int i = 1; i <= y; i++)
	{
		product *= x;
	}
	return product;
}
}

namespace ParticleTracking
{

// Constructor

CollimatorWakeProcess::CollimatorWakeProcess(int modes, int prio, size_t nb, double ns) :
	WakeFieldProcess(prio, nb, ns)
{
	nmodes = modes;
	Cm = new double*[modes + 1];
	Sm = new double*[modes + 1];
	wake_ct = new double*[modes + 1];
	wake_st = new double*[modes + 1];
	wake_cl = new double*[modes + 1];
	wake_sl = new double*[modes + 1];
	for(int i = 0; i <= nmodes; i++)
	{
		Cm[i] = new double[1000];
		Sm[i] = new double[1000];
		wake_ct[i] = new double[1000];
		wake_st[i] = new double[1000];
		wake_cl[i] = new double[1000];
		wake_sl[i] = new double[1000];
	}
}

// Destructor

CollimatorWakeProcess::~CollimatorWakeProcess()
{
	for(int i = 0; i < nmodes; i++)
	{
		delete[] Cm[i];
		delete[] Sm[i];
		delete[] wake_ct[i];
		delete[] wake_st[i];
		delete[] wake_cl[i];
		delete[] wake_sl[i];
	}

	delete[] Cm;
	delete[] Sm;
	delete[] wake_ct;
	delete[] wake_st;
	delete[] wake_cl;
	delete[] wake_sl;
}

// Calculates the moments Cm for each slice
double CollimatorWakeProcess::CalculateCm(int mode, int slice)
{
	double x = 0;
	for(ParticleBunch::iterator p = bunchSlices[slice]; p != bunchSlices[slice + 1]; p++)
	{
		double r = sqrt(powd(p->x(), 2) + powd(p->y(), 2));
		double theta = atan2(p->y(), p->x());
		x += powd(r, mode) * cos(mode * theta);
	}
	return x;
}

// Calculates  the moments Sm for each slice
double CollimatorWakeProcess::CalculateSm(int mode, int slice)
{
	double x = 0;
	for(ParticleBunch::iterator p = bunchSlices[slice]; p != bunchSlices[slice + 1]; p++)
	{
		double r = sqrt(powd(p->x(), 2) + powd(p->y(), 2));
		double theta = atan2(p->y(), p->x());
		x += powd(r, mode) * sin(mode * theta);
	}
	return x;
}

// Calculate the transverse wake with modes
void CollimatorWakeProcess::CalculateWakeT(double dz, int currmode)
{
	vector<double> w(nbins);
	for(size_t slice = 0; slice < nbins; slice++)
	{
		w[slice] = collimator_wake->Wtrans(slice * dz, currmode);
	}

	for(size_t slice = 0; slice < nbins; slice++)
	{
		int i = slice;
		wake_ct[currmode][i] = 0;
		wake_st[currmode][i] = 0;

		for(size_t j = slice; j < bunchSlices.size() - 1; j++)
		{
			wake_ct[currmode][i] += w[j - i] * Cm[currmode][j];
			wake_st[currmode][i] += w[j - i] * Sm[currmode][j];
		}
	}
}

// This function calculates the longitudinal wake with modes
void CollimatorWakeProcess::CalculateWakeL(double dz, int currmode)
{
	vector<double> w(nbins);
	for(size_t slice = 0; slice < nbins; slice++)
	{
		w[slice] = collimator_wake->Wlong(slice * dz, currmode);
	}

	for(size_t i = 0; i < nbins; i++)
	{
		wake_cl[currmode][i] = 0;
		wake_sl[currmode][i] = 0;
		for(size_t j = i; j < bunchSlices.size() - 1; j++)
		{
			wake_cl[currmode][i] += w[j - i] * Cm[currmode][j];
			wake_sl[currmode][i] += w[j - i] * Sm[currmode][j];
		}
	}
}

void CollimatorWakeProcess::ApplyWakefield(double ds) //  int nmodes)
{
	collimator_wake = (CollimatorWakePotentials *) currentWake;
	for(int m = 1; m <= nmodes; m++)
	{
		for(size_t n = 0; n < nbins; n++)
		{
			Cm[m][n] = CalculateCm(m, n);
			Sm[m][n] = CalculateSm(m, n);
		}
	}

	double wake_x, wake_y, wake_z;
	double macrocharge = currentBunch->GetTotalCharge() / currentBunch->size();
	double a0 = macrocharge * ElectronCharge * Volt;
	a0 /= 4 * pi * FreeSpacePermittivity;
	double p0 = currentBunch->GetReferenceMomentum();

	if(recalc)
	{
		Init();
	}
	double bload = 0;

#define WAKE_GRADIENT(wake) ((wake[currmode][nslice + 1] - wake[currmode][nslice]) / dz)

	for(int currmode = 1; currmode <= nmodes; currmode++)
	{
		CalculateWakeT(dz, currmode);
		CalculateWakeL(dz, currmode);
		double z = zmin;
		int iparticle = 0;

		for(size_t nslice = 0; nslice < nbins; nslice++)
		{
			double g_ct = WAKE_GRADIENT(wake_ct);
			double g_st = WAKE_GRADIENT(wake_st);
			double g_cl = WAKE_GRADIENT(wake_cl);
			double g_sl = WAKE_GRADIENT(wake_sl);
			g_ct = g_st = g_cl = g_sl = 0;
			int number_particles = 0;
			for(ParticleBunch::iterator p = bunchSlices[nslice]; p != bunchSlices[nslice + 1]; p++)
			{
				number_particles++;
				double r = sqrt(powd(p->x(), 2) + powd(p->y(), 2));
				double theta = atan2(p->y(), p->x());
				double zz = p->ct() - z;
				double wxc = cos((currmode - 1) * theta) * (wake_ct[currmode][nslice] + g_ct * zz);
				double wxs = sin((currmode - 1) * theta) * (wake_st[currmode][nslice] + g_st * zz);
				double wys = cos((currmode - 1) * theta) * (wake_st[currmode][nslice] + g_st * zz);
				double wyc = sin((currmode - 1) * theta) * (wake_ct[currmode][nslice] + g_ct * zz);
				wake_x = currmode * powd(r, currmode - 1) * (wxc + wxs);
				wake_y = currmode * powd(r, currmode - 1) * (wys - wyc);
				wake_x *= a0;
				wake_y *= a0;
				double wzc = cos(currmode * theta) * (wake_cl[currmode][nslice] + g_cl * zz);
				double wzs = sin(currmode * theta) * (wake_sl[currmode][nslice] + g_sl * zz);
				wake_z = powd(r, currmode) * (wzc - wzs);
				wake_z *= a0;
				double ddp = -wake_z / p0;
				p->dp() += ddp;
				bload += ddp;
				double dxp = inc_tw ? wake_x / p0 : 0;
				double dyp = inc_tw ? wake_y / p0 : 0;
				p->xp() = (p->xp() + dxp) / (1 + ddp);
				double oldpy = p->yp();
				p->yp() = (p->yp() + dyp) / (1 + ddp);
			}

			z += dz;
		}
	}
}
} //end namespace ParticleTracking
