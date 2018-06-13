/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "WakePotentials.h"

#include "CollimatorWakeProcess.h"
#include "CollimatorWakePotentials.h"
#include "CollimatorTable.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

//------------------------------------------------------------------------------------------------------------
//      the resistive wake potentials  (in MKS system)
//------------------------------------------------------------------------------------------------------------

class ResistivePotential: public CollimatorWakePotentials
{
public:
	double sigma, b, leng, scale, step;
	int ncoeff;
	double lo, hi, * coeff;
	collimatortable** Transverse;
	collimatortable** Longitudinal;

	ResistivePotential(int m, double ss, double bb, double l, string filename, double tau = 0) :
		CollimatorWakePotentials(m, 0., 0.), sigma(ss), b(bb), leng(l)
	{
		double Z0 = 377;
		scale = pow(2 * b * b / (Z0 * sigma), 1. / 3.);

		double xi = pow(scale / b, 2);
		double Gamma = SpeedOfLight * tau / scale;
		Transverse = new collimatortable*[m + 1];
		Longitudinal = new collimatortable*[m + 1];
		for(int mode = 0; mode <= m; mode++)
		{
			stringstream ss;
			ss << filename << "L" << mode << ".txt";
			Longitudinal[mode] = new collimatortable(ss.str().c_str(), Gamma, xi);
			if(mode > 0)
			{
				stringstream ss;
				ss << filename << "T2m" << mode << ".txt";
				Transverse[mode] = new collimatortable(ss.str().c_str(), Gamma, xi);
			}
		}
	} // end of constructor

	double Wlong(double z) const
	{
		return 1;
	}
	double Wtrans(double z) const
	{
		return 1;
	}
	double Wtrans(double z, int m) const
	{
		if(z < 0)
		{
			return 0;
		}
		double s = z / scale;
		double fourpieps = 1; // pu in later
		double Chao = -2 * (1 / (sqrt(2 * pi))) * sqrt(scale / z);
		double E = Transverse[m]->inrange(s) ? Transverse[m]->interpolate(s) : Chao;
		E = scale * leng * E; // minus sign should have been in Mathematica
		E = E / fourpieps; // minus sign should have been in Mathematica
		E = E / pow(b, 2 * m + 2);
		return E;
	}

	double Wlong(double z, int m) const
	{
		return z > 0 ? 1 : 0;
	}
}; //End of ResistivePotential Class

class ResistiveWakePotentials: public CollimatorWakePotentials
{
public:

	double* coeff;
	double rad, sigma, length;

	ResistiveWakePotentials(int m, double r, double s, double l) :
		CollimatorWakePotentials(m, r, s), rad(r), sigma(s), length(l)
	{
		cout << "Making new ResistiveWakePotentials with length: " << length << endl;
		coeff = new double[m + 1];

		int delta;
		if(m == 0)
		{
			delta = 1;
		}

		else
		{
			delta = 0;
		}

		for(int i = 0; i < (m + 1); i++)
		{
			coeff[i] = 1 / (pi * pow(rad, 2 * i + 1) * (1 + delta));
		}
	}

	double Wlong(double z) const
	{
		return 1;
	}
	double Wtrans(double z) const
	{
		return 1;
	}
	double Wtrans(double z, int m) const
	{
		cout << " call for " << m << endl;
		return z > 0 ? coeff[m] * sqrt(376.74 / (pi * sigma)) * length / sqrt(z) : 0;
	}
	double Wlong(double z, int m) const
	{
		return z > 0 ? -coeff[m] * sqrt(1 / sigma * 376.6) * sqrt(z) * length : 0;
	}

};
