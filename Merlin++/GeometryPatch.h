/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef GeometryPatch_h
#define GeometryPatch_h 1

#include "merlin_config.h"
#include <utility>

#include "Transformable.h"
#include "AcceleratorGeometry.h"

/**
 * A geometry patch, representing an arbitrary
 * transformation of the accelerator geometry.
 * GeometryPatch objects have zero extents.
 *
 * Note: this class was developed primarily to
 *       support MAD-like SROT elements.
 */

class GeometryPatch: public AcceleratorGeometry, public Transformable
{
public:

	virtual Transform3D GetGeometryTransform(double s0, double s) const
	{
		if(s0 != 0 || s != 0)
		{
			throw BeyondExtent();
		}
		return local_T ? *local_T : Transform3D();
	}

	virtual Transform3D GetGeometryTransform(BoundaryPlane p) const
	{
		return (p == entrance && local_T) ? *local_T : Transform3D();
	}

	virtual Transform3D GetTotalGeometryTransform() const
	{
		return local_T ? *local_T : Transform3D();
	}

	/**
	 *	Returns the local extent of this geometry.
	 *	@return Local extent of the geometry
	 */
	virtual Extent GetGeometryExtent() const
	{
		return Extent(0, 0);
	}

	/**
	 *	Returns the total arc-length of the geometry.
	 *	@return Total arc-length of the geometry
	 */
	virtual double GetGeometryLength() const
	{
		return 0.0;
	}
};

#endif
