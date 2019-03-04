/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Transformable_h
#define Transformable_h 1

#include "merlin_config.h"
#include "Transform3D.h"

/**
 * A mixin class which can be used to add translation and rotation
 * methods to other (derived) classes.
 */
class Transformable
{
public:

	/**
	 * Construction
	 */
	Transformable();
	Transformable(const Transformable& rhs);

	/**
	 * Destruction
	 */
	virtual ~Transformable();

	/**
	 *	Returns true if this object has been locally transformed.
	 */
	bool IsTransformed() const;

	/**
	 * Translate by the vector X
	 */
	void Translate(const Vector3D& X);

	/**
	 * Translate by the vector (x,y z)
	 */
	void Translate(double x, double y, double z);

	/**
	 *	Translates the frame along the current x-axis by dx.
	 */
	void TranslateX(double dx);

	/**
	 *	Translates the frame along the current y-axis by dy.
	 */
	void TranslateY(double dy);

	/**
	 *	Translates the frame along the current z-axis by dz.
	 */
	void TranslateZ(double dz);

	/**
	 *	Rotates the frame about the current x-axis by angle.
	 */
	void RotateX(double angle);

	/**
	 *	Rotates the frame about the current y-axis by angle.
	 */
	void RotateY(double angle);

	/**
	 *	Rotates the frame about the current z-axis by angle.
	 */
	void RotateZ(double angle);

	/**
	 *	Clear the transformation.
	 */
	void ClearTransform();

	/**
	 * Virtual function to allow derived classes to act on
	 * change of state.
	 */
	virtual void Invalidate() const
	{
		/* do nothing */
	}

	const Transform3D* GetTransformation() const
	{
		return local_T;
	}

protected:

	/**
	 *	The local transformation
	 */
	Transform3D* local_T;

	/**
	 * Copy
	 */
	Transformable& operator=(const Transformable& rhs)
	{
		if(local_T)
		{
			delete local_T;
		}
		local_T = rhs.local_T ? new Transform3D(*rhs.local_T) : nullptr;
		return *this;
	}
};

/**
 * Class Transformable
 */

inline Transformable::Transformable() :
	local_T(nullptr)
{
}

inline Transformable::Transformable(const Transformable& rhs) :
	local_T(nullptr)
{
	local_T = rhs.local_T ? new Transform3D(*rhs.local_T) : nullptr;
}

inline bool Transformable::IsTransformed() const
{
	return local_T != nullptr && !local_T->isIdentity();
}

inline void Transformable::Translate(const Vector3D& X)
{
	Translate(X.x, X.y, X.z);
}

inline void Transformable::TranslateX(double dx)
{
	Translate(dx, 0, 0);
}

inline void Transformable::TranslateY(double dy)
{
	Translate(0, dy, 0);
}

inline void Transformable::TranslateZ(double dz)
{
	Translate(0, 0, dz);
}

#endif
