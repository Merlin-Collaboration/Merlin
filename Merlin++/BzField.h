/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BzField_h
#define BzField_h 1

#include "merlin_config.h"
#include "EMField.h"

/**
 *	Represents a constant magnetic field in along the local
 *	z-axis (a solenoidal field.)
 */

class BzField: public EMField
{
public:

	explicit BzField(double B);

	double GetStrength() const;

	/**
	 *	Returns the magnetic field at the point x and time t.
	 *	@return Magnetic field at point x, time t
	 */
	virtual Vector3D GetBFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Returns the electric field at the point x and time t
	 *	@return Electric field at point x, time t
	 */
	virtual Vector3D GetEFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Sets the strength of the field in Tesla.
	 */
	void SetStrength(double B);

private:

	double Bz;
};

inline BzField::BzField(double B) :
	Bz(B)
{
}

inline double BzField::GetStrength() const
{
	return Bz;
}

inline void BzField::SetStrength(double B)
{
	Bz = B;
}

#endif
