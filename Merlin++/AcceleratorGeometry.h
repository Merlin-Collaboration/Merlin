/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorGeometry_h
#define AcceleratorGeometry_h 1

#include "merlin_config.h"
#include <utility>

#include "Transform3D.h"

/**
 * Represents a frame of reference for a section of
 * accelerator lattice. An AcceleratorGeometry can be
 * considered a type of three-dimensional space line
 * (R(s)), which is characterised by a single scalar s, the
 * distance along the space line from the origin. Each
 * AcceleratorGeometry has a specific length which bounds
 * the allowed values of s (with respect to the local
 * geometry origin). At each position s on the geometry,
 * a local rectangular coordinate frame can be uniquely
 * defined, with its origin at the point s, and its z-axis
 * tangential to the geometry at s. The orientation of the
 * local x- and y-axes are also uniquely determined by the
 * sum of any rotations applied going from the origin to s.
 *
 * The primary responsibility for an AcceleratorGeometry
 * object is to supply transformations between coordinate
 * frames defined on that geometry.
 */
class AcceleratorGeometry
{
public:

	/**
	 * Exception thrown indicating that an s-distance was
	 * outside of the current geometry extent.
	 */
	class BeyondExtent
	{
	};

	typedef std::pair<double, double> Extent;

	/**
	 * A BoundaryPlane is the X-Y plane (z=0) of the coordinate
	 * frame defined at the entrance (start) or exit (end) of
	 * the Geometry.
	 */
	typedef enum
	{
		entrance,
		exit

	} BoundaryPlane;

	/**
	 * Virtual destructor.
	 */
	virtual ~AcceleratorGeometry();

	/**
	 * Return the three-dimensional transformation from the
	 * frame at s0 to the frame at s. s and s0 are in the
	 * geometry's s-frame, and must be within the geometry
	 * extents.
	 * @param[in] s0 The location at which the transform should be evaluated from.
	 * @param[in] s The location at which the transform should be evaluated to.
	 * @exception Throws a BeyondExtent exception if the requested s values are outside the geometry extent.
	 * @return The 3D transformation from the entrance to the exit of this geometry.
	 */
	virtual Transform3D GetGeometryTransform(double s0, double s) const = 0;

	/**
	 * Return the three-dimensional transformation from the
	 * local origin to the frame at s. s is in the geometry's
	 * s-frame, and must be within the geometry extents.
	 * @param[in] s The location at which the transform should be evaluated to.
	 * @exception Throws a BeyondExtent exception if the requested s value is outside the geometry extent.
	 * @return The 3D transformation from the entrance to the exit of this geometry.
	 */
	Transform3D GetGeometryTransform(double s) const;

	/**
	 * Returns the transformation from the geometry origin to
	 * the specified boundary plane.
	 * @param[in] p The chosen BoundaryPlane
	 * @return The 3D transformation from the geometry origin to a specified boundary plane (entrance or exit).
	 */
	virtual Transform3D GetGeometryTransform(BoundaryPlane p) const = 0;

	/**
	 * Returns the transformation from the entrance plane frame
	 * to the exit plane frame.
	 * @return The 3D transformation from the entrance to the exit of this geometry.
	 */
	virtual Transform3D GetTotalGeometryTransform() const;

	/**
	 * Returns the local extent of this geometry.
	 * @return A pair giving the extent of this geometry.
	 */
	virtual Extent GetGeometryExtent() const = 0;

	/**
	 * Returns the total arc-length of the geometry.
	 * @return A double giving the total length of this geometry.
	 */
	virtual double GetGeometryLength() const = 0;
};

inline AcceleratorGeometry::~AcceleratorGeometry()
{
}

inline Transform3D AcceleratorGeometry::GetGeometryTransform(double s) const
{
	return GetGeometryTransform(0, s);
}

inline double ToLength(const AcceleratorGeometry::Extent& extent)
{
	return extent.second - extent.first;
}

#endif
