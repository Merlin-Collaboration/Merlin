/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MagnetMover_h
#define MagnetMover_h 1

#include "merlin_config.h"
#include "SequenceFrame.h"
#include "Transform2D.h"

/**
 *	Represents a mechanical magnet mover, or a remotely
 *	translatable stage.  A magnet mover can effectively
 *	apply a 2-D transformation with respect to its local
 *	reference frame. The transformation is characterised by
 *	a rotation about the local z-axis (roll) followed by a
 *	translation in the (x-y) plane.
 *
 *	The local coordinate frame can be rotated and translated
 *	using the methods inherited from LatticeFrame. The
 *	additional 2-D mover transformations are in addition to
 *	any local (3-D) frame modification. In this way, magnet
 *	mover alignment errors can be simulated.
 */

class MagnetMover: public SequenceFrame
{
public:

	/**
	 *	Constructor taking the name of the mover, and the Lattice
	 *	Frame which the mover adjusts.
	 */
	explicit MagnetMover(const string& id);

	/**
	 *	Returns the horizontal displacement in meters.
	 *
	 *	@return Horizontal displacement (m)
	 */
	double GetX() const;

	/**
	 *	Returns the vertical displacement in meters.
	 *
	 *	@return Vertical displacement (m)
	 */
	double GetY() const;

	/**
	 *	Returns the roll angle about the local entrance z-axis
	 *	in radians.
	 *
	 *	@return Roll angle (radians)
	 */
	double GetRoll() const;

	/**
	 *	Sets the horizontal offset in meters.
	 */
	void SetX(double x);

	/**
	 *	Sets the vertical offset in meters.
	 */
	void SetY(double y);

	/**
	 *	Sets the roll angle in radians.
	 */
	void SetRoll(double roll);

	/**
	 *	Resets the mover to zero.
	 */
	void Reset();

	/**
	 *	Returns MagnetMover
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor
	 */
	virtual ModelElement* Copy() const;

	/**
	 *	Returns the total frame transformation, i.e. the
	 *	combined transformation of the local frame transforms
	 *	and the x, y and roll  mover adjustments.
	 */
	virtual Transform3D GetLocalFrameTransform() const;

private:

	Transform2D t;
};

inline MagnetMover::MagnetMover(const string& id) :
	SequenceFrame(id)
{
}

inline double MagnetMover::GetX() const
{
	return t.translation().x;
}

inline double MagnetMover::GetY() const
{
	return t.translation().y;
}

inline double MagnetMover::GetRoll() const
{
	return t.rotationAngle();
}

inline void MagnetMover::SetX(double x)
{
	t.setTranslationX(x);
}

inline void MagnetMover::SetY(double y)
{
	t.setTranslationY(y);
}

inline void MagnetMover::SetRoll(double roll)
{
	t.setRotation(roll);
}

inline void MagnetMover::Reset()
{
	t = Transform2D();
}

#endif
