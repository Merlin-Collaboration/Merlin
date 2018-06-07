/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorSupport_h
#define AcceleratorSupport_h 1

#include "merlin_config.h"
#include <vector>
#include "Space2D.h"
#include "Space3D.h"

/**
 *	Represents a single support structure. Support can be
 *	translated in either x, y or z directions. Support is
 *	primarily intended for ground motion application.
 */

/*
 *	An array of AcceleratorSupport pointers.
 */
class AcceleratorSupport;
typedef std::vector<AcceleratorSupport*> AcceleratorSupportList;

class AcceleratorSupport
{
public:

	AcceleratorSupport();

	/**
	 *	Sets the location of the support in the global
	 *	coordinate frame.
	 *	@param[in] x The x position in the accelerator plane
	 *	@param[in] z The z position in the accelerator plane
	 *	@param[in] s The arc position of the support
	 */
	void SetPosition(double s, double x, double z);

	/*
	 *	Returns the arc position.
	 *	@return The arc position
	 */
	double GetArcPosition() const;

	/**
	 *	Returns the location of the support in the accelerator
	 *	plane (x,z). Note that  Point2D::y here refers to the
	 *	z-coordinate.
	 *	@return The location of the support in the accelerator plane
	 */
	Point2D GetLocation() const;

	/*
	 *	Returns the linear distance from this support to another Support.
	 *	@return The linear distance from this support to another support
	 */
	double DistanceTo(const AcceleratorSupport& aSupport) const;

	/*
	 *	Set the support offset to (x,y,z).
	 */
	void SetOffset(double x, double y, double z);

	/*
	 *	Set the support offset to X
	 */
	void SetOffset(const Vector3D& X);

	/*
	 *	Returns the current offset.
	 *	@return The current offset
	 */
	const Vector3D& GetOffset() const;

	/*
	 *	Increment the current offset by (dx,dy,dz),
	 */
	const Vector3D& IncrementOffset(double dx, double dy, double dz);

	/*
	 *	Increment the current offset by dX.
	 */
	const Vector3D& IncrementOffset(const Vector3D& dX);

	/*
	 *	Reset the offset to (0,0,0).
	 */
	void Reset();

private:

	Vector3D offset;
	bool modified;
	Point2D pos;
	double s_pos;
	friend class SupportStructure;
};

inline AcceleratorSupport::AcceleratorSupport() :
	offset(0, 0, 0), modified(false), pos(0, 0), s_pos(0)
{
}

inline void AcceleratorSupport::SetPosition(double s, double x, double z)
{
	pos.x = x;
	pos.y = z;
	s_pos = s;
}

inline double AcceleratorSupport::GetArcPosition() const
{
	return s_pos;
}

inline Point2D AcceleratorSupport::GetLocation() const
{
	return pos;
}

inline void AcceleratorSupport::SetOffset(double x, double y, double z)
{
	SetOffset(Vector3D(x, y, z));
}

inline void AcceleratorSupport::SetOffset(const Vector3D& X)
{
	offset = X;
	modified = true;
}

inline const Vector3D& AcceleratorSupport::GetOffset() const
{
	return offset;
}

inline const Vector3D& AcceleratorSupport::IncrementOffset(double dx, double dy, double dz)
{
	return IncrementOffset(Vector3D(dx, dy, dz));
}

inline const Vector3D& AcceleratorSupport::IncrementOffset(const Vector3D& dX)
{
	modified = true;
	return offset += dX;
}

inline void AcceleratorSupport::Reset()
{
	modified = true;
	offset = Vector3D(0, 0, 0);
}

#endif
