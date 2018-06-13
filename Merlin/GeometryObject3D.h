/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef GeometryObject3D_h
#define GeometryObject3D_h 1

#include "merlin_config.h"
#include "Transform3D.h"

/**
 *	A mixin class which represents an geometric entity which
 *	can be rotated and translated with respect to some
 *	(implicit) reference frame.
 */
class GeometryObject3D
{
public:

	/**
	 *	Virtual constructor.
	 */
	virtual ~GeometryObject3D();

	/**
	 *	Translate the object by the relative vector (dx,dy,dz).
	 */
	void Translate(double dx, double dy, double dz);

	/**
	 *	Rotates the object about the current x-axis by angle.
	 */
	void RotateX(double angle);

	/**
	 *	Rotates the object about the current y-axis by angle.
	 */
	void RotateY(double angle);

	/**
	 *	Rotates the object about the current z-axis by angle.
	 */
	void RotateZ(double angle);

	/**
	 *	Transform the object (with respect to the current axes)
	 *	by the transformation t.
	 */
	void ApplyTransformation(const Transform3D& t1);

	/**
	 *	Clear all transformations.
	 */
	void ClearTransformation();

	/**
	 *	Set the absolute transformation for this object.
	 */
	void SetTransformation(const Transform3D& t1);

	/**
	 *	Return the current transformation.
	 */
	const Transform3D& GetTransformation() const;

	/**
	 *	Returns true if this object has been transformed.
	 *	@retval true Object has been transformed
	 *	@retval false Object hasn't been transformed
	 */
	bool IsTransformed() const;

protected:

	/**
	 *	Protected constructor taking the initial transformation.
	 */
	explicit GeometryObject3D(const Transform3D& t0 = Transform3D());

	/**
	 *	Virtual function used to notify derived classes that the
	 *	state of the transformation has changed. Default action
	 *	does nothing.
	 */
	virtual void HasChanged() const;

private:

	Transform3D t;
};

inline GeometryObject3D::GeometryObject3D(const Transform3D& t0) :
	t(t0)
{
}

inline GeometryObject3D::~GeometryObject3D()
{
}

inline void GeometryObject3D::Translate(double dx, double dy, double dz)
{
	t *= Transform3D::translation(dx, dy, dz);
	HasChanged();
}

inline void GeometryObject3D::RotateX(double angle)
{
	t *= Transform3D::rotationX(angle);
	HasChanged();
}

inline void GeometryObject3D::RotateY(double angle)
{
	t *= Transform3D::rotationY(angle);
	HasChanged();
}

inline void GeometryObject3D::RotateZ(double angle)
{
	t *= Transform3D::rotationZ(angle);
	HasChanged();
}

inline void GeometryObject3D::ApplyTransformation(const Transform3D& t1)
{
	t *= t1;
	HasChanged();
}

inline void GeometryObject3D::ClearTransformation()
{
	t = Transform3D();
	HasChanged();
}

inline void GeometryObject3D::SetTransformation(const Transform3D& t1)
{
	t = t1;
	HasChanged();
}

inline const Transform3D& GeometryObject3D::GetTransformation() const
{
	return t;
}

inline bool GeometryObject3D::IsTransformed() const
{
	return !t.isIdentity();
}

inline void GeometryObject3D::HasChanged() const
{
}

#endif
