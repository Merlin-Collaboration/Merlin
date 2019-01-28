/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SWRFfield_h
#define SWRFfield_h 1

#include "merlin_config.h"
#include <cmath>
#include "RFAcceleratingField.h"
#include "PhysicalConstants.h"

/**
 *	A RF standing wave  field as found in standing wave
 *	accelerating structures. The field (0,0,Ez) is defined as
 */
//	        Ez(z,t)=E0 cos(k*z)cos(w*t+phi)
/*
 *   \f[
 *           E_z(z,t)=E_0 \cos(kz)\cos(\omega t+\phi)
 *   \f]
 */

class SWRFfield: public RFAcceleratingField
{
public:

	/**
	 *	Constructor taking the frequency, peak electric field
	 *	and the phase of the RF.
	 */
	SWRFfield(double f, double Epk, double phase = 0);

	/**
	 *	Returns the magnetic field at the point x and time t.
	 *	@return Magnetic field at (x,t)
	 */
	virtual Vector3D GetBFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Returns the electric field at the point x and time t
	 *	@return Electric field at (x,t)
	 */
	virtual Vector3D GetEFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Returns the force due to this field on a particle of
	 *	charge q with position x and velocity v at time t.
	 *
	 *	@return Force acting on charged particle
	 */
	virtual Vector3D GetForceAt(const Point3D& x, const Vector3D& v, double q, double t = 0) const;

	/**
	 *	Calculate the Ez component.
	 */
	virtual double Ez(double z, double t) const;
};

inline SWRFfield::SWRFfield(double f, double Epk, double phase) :
	RFAcceleratingField(Epk, f, phase)
{
}

inline double SWRFfield::Ez(double z, double t) const
{
	using PhysicalConstants::SpeedOfLight;
	return E0 * cos(w * z / SpeedOfLight) * cos(w * t + phi);
}

#endif
