/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BeamData_h
#define BeamData_h 1

#include "merlin_config.h"
#include <cmath>

/**
 *	Data structure for defining the 6-D beam phase space
 *	(first- and second-order moments).
 */

class BeamData
{
public:

	BeamData() :
		beta_x(0), beta_y(0), alpha_x(0), alpha_y(0), emit_x(0), emit_y(0), sig_dp(0), sig_z(0), x0(0), xp0(0), y0(0),
		yp0(0), ct0(0), p0(0), c_xy(0), c_xyp(0), c_xpy(0), c_xpyp(0), Dx(0), Dxp(0), Dy(0), Dyp(0), charge(1)
	{
	}

	/**
	 *	Checks consistency of definition. Returns true if *this
	 *	is a valid beam, otherwise false.
	 *	@retval true If beam is valid
	 *	@retval false If beam is invalid
	 */
	bool ok() const
	{
		return emit_x > 0 && emit_y > 0 && beta_x > 0 && beta_y > 0;
	}

	// TWISS beam parameters
	double beta_x;
	double beta_y;
	double alpha_x;
	double alpha_y;
	double gamma_x() const
	{
		return (1 + alpha_x * alpha_x) / beta_x;
	}
	double gamma_y() const
	{
		return (1 + alpha_y * alpha_y) / beta_y;
	}
	double emit_x;
	double emit_y;

	/**
	 *	Relative energy spread of beam.
	 */
	double sig_dp;

	/**
	 *	Beam length.
	 */
	double sig_z;

	//	Beam centroid.
	double x0;
	double xp0;
	double y0;
	double yp0;
	double ct0;

	/**
	 *	Beam energy (momentum).
	 */
	double p0;

	//	X-Y coupling
	double c_xy;
	double c_xyp;
	double c_xpy;
	double c_xpyp;

	//	Dispersion
	double Dx;
	double Dxp;
	double Dy;
	double Dyp;

	/**
	 *	The charge of the particles in the beam.
	 *  <0 for electrons, >0 for positrons.
	 */
	double charge;

	double sigma_x() const
	{
		return sqrt(emit_x * beta_x);
	}
	double sigma_y() const
	{
		return sqrt(emit_y * beta_y);
	}
	double sigma_xp() const
	{
		return sqrt(emit_x * gamma_x());
	}
	double sigma_yp() const
	{
		return sqrt(emit_y * gamma_y());
	}
};

#endif
