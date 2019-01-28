/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CenteredGeometry_h
#define CenteredGeometry_h 1

#include "merlin_config.h"
#include "AcceleratorGeometry.h"

/**
 *	Base class for accelerator geometry types whose origin
 *	is at the center of the their arc length. Hence their
 *	extent extends from -length/2 to +length/2.
 */

class CenteredGeometry: public AcceleratorGeometry
{
public:

	/**
	 *	Constructor taking the arc length of the geometry.
	 */
	explicit CenteredGeometry(double l);

	/**
	 *	Returns the total arc-length of the geometry.
	 *	@return Total arc-length
	 */
	virtual double GetGeometryLength() const;

	/**
	 *	Returns the local extent of this geometry.
	 *	@return Local extent of this geometry
	 */
	virtual AcceleratorGeometry::Extent GetGeometryExtent() const;

protected:

	/**
	 *	Used to check if (s1s2) is within the geometry bounds.
	 */
	void CheckBounds(double s1, double s2) const;

	/**
	 *	Used to check if s1 is within the geometry bounds.
	 */
	void CheckBounds(double s) const;

	double len;
};

inline CenteredGeometry::CenteredGeometry(double l) :
	AcceleratorGeometry(), len(l)
{
}

inline double CenteredGeometry::GetGeometryLength() const
{
	return len;
}

inline AcceleratorGeometry::Extent CenteredGeometry::GetGeometryExtent() const
{
	double s = len / 2;
	return Extent(-s, s);
}

inline void CenteredGeometry::CheckBounds(double s1, double s2) const
{
	double hl = len / 2;
	if(fabs(s1) > hl || fabs(s2) > hl)
	{
		throw BeyondExtent();
	}
}

inline void CenteredGeometry::CheckBounds(double s) const
{
	if(fabs(s) > len / 2)
	{
		throw BeyondExtent();
	}
}

#endif
