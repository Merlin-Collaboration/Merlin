/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RectangularGeometry_h
#define RectangularGeometry_h 1

#include "merlin_config.h"
#include "CenteredGeometry.h"

/**
 *	Represents a straight line segment. Transformations from
 *	points on a RectangularGeometry are pure translations
 *	along the z-axis.
 */

class RectangularGeometry: public CenteredGeometry
{
public:

	/**
	 *	Constructor taking the length of the rectangular
	 *	geometry (total z extent).
	 */
	RectangularGeometry(double l);

	/**
	 *	Returns a translation along the z-axis of (s-s0).
	 *	@return Translation along z of \f$ s-s_0 \f$
	 */
	virtual Transform3D GetGeometryTransform(double s0, double s) const;

	/**
	 *	Returns a translation along the z-axis of either +l/2 or
	 *	-l/2 for the entrance and exit boundary planes
	 *	respectively.
	 *
	 *	@return Translation along z of \f$ \pm\frac{1}{2} \f$ corresponding to
	 *	the entrance/exit boundary planes respectively
	 */
	virtual Transform3D GetGeometryTransform(BoundaryPlane p) const;

	/**
	 *	Returns a translation along the z-axis of +l.
	 *	@return Translation along z of \f$ +1 \f$
	 */
	virtual Transform3D GetTotalGeometryTransform() const;
};

inline RectangularGeometry::RectangularGeometry(double l) :
	CenteredGeometry(l)
{
}

#endif
