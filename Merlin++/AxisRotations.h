/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AxisRotations_h
#define AxisRotations_h 1

#include "merlin_config.h"
#include <math.h>
#include "Space3D.h"
#include "RotationType.h"
#include "Rot3Drep.h"

class GeneralRotation;

/**
 *	Base class for pure rotations about a specific axis.
 */

class PureAxisRotation: public Rot3Drep
{
public:

	/**
	 *	Constructor taking the axis rotation angle in radians.
	 *	@param[in] phi Axis rotation angle (radians)
	 */
	PureAxisRotation(double phi);

	/**
	 *	Return the angle of rotation in radians.
	 *	@return Angle of rotation (radians)
	 */
	double angle() const;

	/**
	 *	Return the sine of the rotation angle.
	 *	@return Sine of the rotation angle
	 */
	double sine() const;

	/**
	 *	Return the cosine of the rotation angle.
	 *	@return Cosine of the rotation angle
	 */
	double cosine() const;

protected:

	/**
	 *	Protected constructor taking the sine and cosine of the
	 *	rotation angle
	 *	@param cosa Cosine of the rotation angle
	 *	@param sina Sine of the rotation angle
	 */
	PureAxisRotation(double cosa, double sina);

private:

	/*
	 *   The single angle is stored as sine and cosine for efficiency reasons when
	 *   performing rotations.
	 */
	double cosphi;
	double sinphi;
};

/**
 * A pure rotation about the x-axis.
 */

class RotationX: public PureAxisRotation
{
public:

	/**
	 *	Constructor taking the rotation angle in radians.
	 *	@param[in] phi Rotation angle (radians)
	 */
	RotationX(double phi);

	/**
	 *	Return an inverted rotation.
	 */
	virtual Rot3Drep* inv() const;

	/**
	 *	Rotate the specified point and return the result.
	 *	@param[in] p A point in 3D space
	 *	@return The rotated point @ref p
	 */
	virtual Point3D rotate(const Point3D& p) const;

	/**
	 *	Rotate the specified vector and return the result.
	 *	@param[in] v A 3D vector
	 *	@return The rotated vector @ref v
	 */
	virtual Vector3D rotate(const Vector3D& v) const;

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
	 *	@return Rotation type
	 */
	virtual RotationType type() const;

	/**
	 *	Returns true.
	 *	@return true
	 */
	virtual bool isXrot() const;

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
	 *	@param[out] m 3x3 rotation matrix
	 */
	virtual RealMatrix& getMatrix(RealMatrix& m) const;

private:

	/**
	 *	Private constructor taking the sine and cosine of the
	 *	rotation angle.
	 *	@param[in] cosa Cosine of rotation angle
	 *	@param[in] sina Sine of rotation angle
	 */
	RotationX(double cosa, double sina);
};

/**
 * A pure rotation about the y-axis.
 */

class RotationY: public PureAxisRotation
{
public:
	/**
	 *	Constructor taking the rotation angle in radians.
	 *	@param[in] phi Rotation angle (radians)
	 */
	RotationY(double phi);

	/**
	 *	Return an inverted rotation.
	 */
	virtual Rot3Drep* inv() const;

	/**
	 *	Rotate the specified point and return the result.
	 *	@param[in] p A point in 3D space
	 *	@return The rotated point @ref p
	 */
	virtual Point3D rotate(const Point3D& p) const;

	/**
	 *	Rotate the specified vector and return the result.
	 *	@param[in] v A 3D vector
	 *	@return The rotated vector @ref v
	 */
	virtual Vector3D rotate(const Vector3D& v) const;

	/**
	 *	Dot this rotation with r.
	 */
	virtual Rot3Drep* dot(const Rot3Drep& r) const;

	/**
	 *	Dot this rotation by a pure x rotation.
	 */
	virtual Rot3Drep* dotBy(const RotationX& rx) const;

	//	Dot this rotation by a pure y rotation.
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
	 *	@return Rotation type
	 */
	virtual RotationType type() const;

	/**
	 *	Returns true.
	 *	@return true
	 */
	virtual bool isYrot() const;

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
	 *	@param[out] m 3x3 rotation matrix
	 */
	virtual RealMatrix& getMatrix(RealMatrix& m) const;

protected:
private:
	/**
	 *	Private constructor taking the sine and cosine of the
	 *	rotation angle.
	 */
	RotationY(double cosa, double sina);

private:
};

/**
 * A pure rotation about the z-axis.
 */
class RotationZ: public PureAxisRotation
{
public:
	/**
	 *	Constructor taking the rotation angle in radians.
	 *	@param[in] phi Axis rotation angle (radians)
	 */
	RotationZ(double phi);

	/**
	 *	Return an inverted rotation.
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
	 *	@return Rotation type
	 */
	virtual RotationType type() const;

	/**
	 *	Returns true.
	 *	@return true
	 */
	virtual bool isZrot() const;

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
	 *	@param[out] m 3x3 rotation matrix
	 */
	virtual RealMatrix& getMatrix(RealMatrix& m) const;

protected:
private:
	/**
	 *	Private constructor taking the sine and cosine of the
	 *	rotation angle.
	 */
	RotationZ(double cosa, double sina);

private:
};

/**
 * Class PureAxisRotation
 */

inline PureAxisRotation::PureAxisRotation(double phi) :
	cosphi(cos(phi)), sinphi(sin(phi))
{
}

inline PureAxisRotation::PureAxisRotation(double cosa, double sina) :
	cosphi(cosa), sinphi(sina)
{
}

inline double PureAxisRotation::sine() const
{
	return sinphi;
}

inline double PureAxisRotation::cosine() const
{
	return cosphi;
}

inline RotationX::RotationX(double phi) :
	PureAxisRotation(phi)
{
}

inline RotationX::RotationX(double cosa, double sina) :
	PureAxisRotation(cosa, sina)
{
}

/**
 * Class RotationY
 */

inline RotationY::RotationY(double phi) :
	PureAxisRotation(phi)
{
}

inline RotationY::RotationY(double cosa, double sina) :
	PureAxisRotation(cosa, sina)
{
}

/**
 * Class RotationZ
 */

inline RotationZ::RotationZ(double phi) :
	PureAxisRotation(phi)
{
}

inline RotationZ::RotationZ(double cosa, double sina) :
	PureAxisRotation(cosa, sina)
{
}

#endif
