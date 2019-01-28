/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CoulombScatter_h
#define CoulombScatter_h 1
#include <cmath>
#include "RandomNG.h"
namespace ParticleTracking
{

/**
 * calculate random small-angle Coulomb scattering
 */
pair<double, double> CoulombScatter(double x, double theta0)
{
	/**
	 * x - material length in meters
	 * theta0 - RMS scattering angle (plane)
	 * See particle data book section 27.3
	 */
	static const double root12 = sqrt(12.0);

	double z1 = RandomNG::normal(0, 1);
	double z2 = RandomNG::normal(0, 1);

	double theta_plane = z2 * theta0;
	double y_plane = z1 * x * theta0 / root12 + x * theta_plane / 2;

	return make_pair(y_plane, theta_plane);
}
}
#endif
