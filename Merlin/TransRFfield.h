/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef TransverseRFfield_h
#define TransverseRFfield_h 1

#include "merlin_config.h"
#include <cmath>
#include "RFAcceleratingField.h"
#include "PhysicalConstants.h"

/**
 *	A RF travelling wave transverse electric field .
 *  The field (Er*cos(theta),Er*sin(theta),0) is defined as
 */
//	        Er(z,t)=E0 cos(k*z-w*t+phi).
/*
 * \f[
 *           E_r(z,t) = E_0 \cos(kz - \omega t + \phi)
 * \f]
 */

class TransverseRFfield: public RFAcceleratingField
{
public:

	/**
	 *	Constructor taking the frequency, peak electric field
	 *	and the phase of the RF. Default orientation is horizontal.
	 */
	TransverseRFfield(double f, double Epk, double phase = 0, double theta = 0);

	/**
	 *	Returns the magnetic field at the point x and time t.
	 *	@return \f$ B(x,t) \f$
	 */
	virtual Vector3D GetBFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Returns the electric field at the point x and time t
	 *	@return \f$ E(x,t) \f$
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
	 *	Calculate the Er component.
	 */
	virtual double Er(double z, double t) const;

	/**
	 * Returns the field orientation for the transverse cavity.
	 * theta=0 represents a purely horizontal field.
	 *
	 * @return The transverse cavity's field orientation
	 * @retval 0 A purely horizontal field
	 */
	double GetFieldOrientation() const
	{
		return theta;
	}
	void SetFieldOrientation(double t)
	{
		theta = t;
	}

	/**
	 * Calculate the Ez component
	 */
	// TODO: TransverseRFfield should not strictly be derived from RFAcceleratingField,
	//       primarily because of this override.
	virtual double Ez(double z, double t) const
	{
		return 0;
	}

private:
	double theta;
};

inline TransverseRFfield::TransverseRFfield(double f, double Epk, double phase, double t) :
	RFAcceleratingField(Epk, f, phase), theta(t)
{
}

inline double TransverseRFfield::Er(double z, double t) const
{
	using PhysicalConstants::SpeedOfLight;
	return E0 * cos(w * (z / SpeedOfLight - t) + phi);
}

#endif
