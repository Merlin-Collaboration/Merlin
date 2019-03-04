/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Rot3Drep_h
#define Rot3Drep_h 1

#include "merlin_config.h"
#include "LinearAlgebra.h"
#include "Space3D.h"
#include "RotationType.h"

class IdentityRotation;
class GeneralRotation;
class RotationZ;
class RotationY;
class RotationX;

/**
 *	The abstract representation (letter class) for Rotation3D
 */

class Rot3Drep
{
public:

	/**
	 *	Virtual destructor.
	 */
	virtual ~Rot3Drep();

	/**
	 *	Rotate the specified point and return the result.
	 */
	virtual Point3D rotate(const Point3D& x) const = 0;

	/**
	 *	Rotate the specified vector and return the result.
	 */
	virtual Vector3D rotate(const Vector3D& v) const = 0;

	/**
	 *	Dot this rotation with r.
	 */
	virtual Rot3Drep* dot(const Rot3Drep& r) const = 0;

	/**
	 *	Dot this rotation by a pure x rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationX& rx) const = 0;

	/**
	 *	Dot this rotation by a pure y rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationY& ry) const = 0;

	/**
	 *	Dot this rotation by a pure z rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationZ& rz) const = 0;

	/**
	 *	Dot this rotation by a general rotation.
	 */
	virtual Rot3Drep* dotBy(const GeneralRotation& r) const = 0;

	/**
	 *	Return an inverted rotation.
	 *	@return Inverted rotation
	 */
	virtual Rot3Drep* inv() const = 0;

	/**
	 *	Return the type of rotation. Can be identity, xrot,
	 *	yrot, zrot or general.
	 *	@retval xrot Rotation type
	 *	@retval yrot Rotation type
	 *	@retval zrot Rotation type
	 *	@retval general Rotation type
	 */
	virtual RotationType type() const = 0;

	/**
	 *	Returns true if this is a null rotation.
	 *	@retval true If a null rotation
	 *	@retval false If not a null rotation
	 */
	virtual bool isIdentity() const;

	/**
	 *	Returns true if a pure rotation about the x-axis.
	 *	@retval true If pure rotation about x
	 *	@retval false If not pure rotation about x
	 */
	virtual bool isXrot() const;

	/**
	 *	Returns true if a pure rotation about the y-axis.
	 *	@retval true If pure rotation about y
	 *	@retval false If not pure rotation about y
	 */
	virtual bool isYrot() const;

	/**
	 *	Returns true if a pure rotation about the z-axis.
	 *	@retval true If pure rotation about z
	 *	@retval false If not pure rotation about z
	 */
	virtual bool isZrot() const;

	/**
	 *	Static factory method. Constructs an identity (null)
	 *	rotation.
	 */
	static Rot3Drep* identity();

	/**
	 *	Static factory method. Constructs a pure x rotation of
	 *	angle radians.
	 */
	static Rot3Drep* rotationX(double angle);

	/**
	 *	Static factory method. Constructs a pure y rotation of
	 *	angle radians.
	 */
	static Rot3Drep* rotationY(double angle);

	/**
	 *	Static factory method. Constructs a pure z rotation of
	 *	angle radians.
	 */
	static Rot3Drep* rotationZ(double angle);

	/**
	 *	Rotation by 180 degrees about the X-axis.
	 */
	virtual Rot3Drep* rotXbyPI() const = 0;

	/**
	 *	Rotation by 180 degrees about the Y-axis.
	 */
	virtual Rot3Drep* rotYbyPI() const = 0;

	/**
	 *	Rotation by 180 degrees about the Z-axis.
	 */
	virtual Rot3Drep* rotZbyPI() const = 0;

	/**
	 *	Returns in m the 3x3 rotation matrix.
	 *	@param[out] m 3x3 rotation matrix
	 */
	virtual RealMatrix& getMatrix(RealMatrix& m) const = 0;

protected:

	/**
	 *	Protected default constructor.
	 */
	Rot3Drep();

private:

	/**
	 *	Reference count for memory management.
	 */
	int refc;

	friend class Rotation3D;
};

#endif
