/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef IdentityRotation_h
#define IdentityRotation_h 1

#include "merlin_config.h"
#include "Space3D.h"
#include "RotationType.h"
#include "Rot3Drep.h"

class GeneralRotation;
class RotationZ;
class RotationY;
class RotationX;

/**
 *	A null or identity rotation.
 */

class IdentityRotation: public Rot3Drep
{
public:

	/**
	 *	Return an inverted rotation.
	 *	@return Inverted rotation
	 */
	virtual Rot3Drep* inv() const;

	/**
	 *	Rotate the specified point and return the result.
	 */
	virtual Point3D rotate(const Point3D& p) const;

	/**
	 *	Rotate the specified vector and return the result.
	 */
	virtual Vector3D rotate(const Vector3D& v) const;

	/**
	 *	Returns true.
	 *	@return true
	 */
	virtual bool isIdentity() const;

	/**
	 *	Dot this rotation with r.
	 */
	virtual Rot3Drep* dot(const Rot3Drep& r) const;

	/**
	 *	Dot this rotation by a pure x rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationX& rx) const;

	/**
	 *	Dot this rotation by a pure y rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationY& ry) const;

	/**
	 *	Dot this rotation by a pure z rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationZ& rz) const;

	/**
	 *	Dot this rotation by a general rotation.
	 */
	virtual Rot3Drep* dotBy(const GeneralRotation& r) const;

	/**
	 *	Return the type of rotation. Can be identity, xrot,
	 *	yrot, zrot or general.
	 *	@return Rotation type (`xrot`, `yrot`, `zrot` or `general`)
	 */
	virtual RotationType type() const;

	/**
	 *	Rotation by 180 degrees about the X-axis.
	 */
	virtual Rot3Drep* rotXbyPI() const;

	/**
	 *	Rotation by 180 degrees about the Y-axis.
	 */
	virtual Rot3Drep* rotYbyPI() const;

	/**
	 *	Rotation by 180 degrees about the Z-axis.
	 */
	virtual Rot3Drep* rotZbyPI() const;

	/**
	 *	Returns in m the 3x3 rotation matrix.
	 *	@return 3x3 rotation matrix
	 */
	virtual RealMatrix& getMatrix(RealMatrix& m) const;
};

#endif
