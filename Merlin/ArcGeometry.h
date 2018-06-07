/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ArcGeometry_h
#define ArcGeometry_h 1

#include "merlin_config.h"
#include "CenteredGeometry.h"

/**
 *	An ArcGeometry represents a constant radius curve in the
 *	local x-z plane. By convention, a positive curvature (h)
 *	defines a curve towards negative x. Transformations
 *	between two points (s1,s2) on an ArcGeometry are
 *	specified by a translation in the x-z plane of
 */
// [(cos(phi)-1)/h, 0, sin(phi)/h], and a rotation about
/*
 *   \f[
 *       [\frac{\cos(\phi-1)}{h}, 0, \frac{\sin(\phi)}{h}]
 *   \f]
 *
 *	and a rotation about the y-axis of
 */
//	-phi (phi = h*(s2-s1)).
/*
 *   /f$ -\phi \quad (\phi=h(s_2-s_1))\f$
 *
 *	An ArcGeometry can have an additional tilt, which is
 *	defined as a rotation about the local z-axis of the
 *	entrance plane by an angle theta, with an additional
 *	rotation about the exit plane z-axis by -theta.
 */

class ArcGeometry: public CenteredGeometry
{
public:

	/**
	 *	Constructor taking the arc length and the curvature
	 *	(1/r) of the geometry.
	 *	@param[in] 1 Arc length
	 *	@param[in] curv Curvature (\f$\frac{1}{r}\f$)
	 */
	ArcGeometry(double l, double curv);

	/**
	 *	Return the curvature of the ArcGeometry.
	 *	@return Curvature of ArcGeometry
	 */

	double GetCurvature() const;

	/**
	 *	Return the total arc angle of the geometry.
	 *	@return Total arc angle of the geometry
	 */
	double GetAngle() const;

	/**
	 *	Returns the arc transform from s0 to s (angle=h*(s-s0)).
	 *	@return Arc transform from s0 to s (\f$\mathrm{angle}=h\times(s-s_0)\f$)
	 */
	virtual Transform3D GetGeometryTransform(double s0, double s) const;

	/**
	 *	Returns the arc transform for either -angle/2 or
	 *	+angle/2 for the entrance and exit planes respectively.
	 *
	 *	@return Arc transform for \f$\frac{-\mathrm{angle}}{2}\f$ or
	 *   \f$\frac{+\mathrm{angle}}{2}\f$ corresponding to the entrance and exit
	 *   planes respectively
	 */
	virtual Transform3D GetGeometryTransform(BoundaryPlane p) const;

	/**
	 *	Returns the arc transformation from the entrance plane
	 *	frame to the exit plane frame.
	 *	@return Arc transformation from the entrance plane frame to the exit
	 *	plane frame
	 */
	virtual Transform3D GetTotalGeometryTransform() const;

	/**
	 *	Sets the curvature (=1/r) of the arc.
	 *	@param[in] curv Curvature (\f$\frac{1}{r}\f$)
	 */
	void SetCurvature(double curv);

	/**
	 *	Sets the tilt of the geometry in radians.
	 */
	void SetTilt(double t);

	/**
	 *	Returns the tilt of the geometry in radians.
	 *	@return Tilt of the geometry (in radians)
	 */
	double GetTilt() const;

private:

	double h;
	double tilt;
};

inline ArcGeometry::ArcGeometry(double l, double curv) :
	CenteredGeometry(l), h(curv), tilt(0)
{
}

inline double ArcGeometry::GetCurvature() const
{
	return h;
}

inline double ArcGeometry::GetAngle() const
{
	return len * h;
}

inline void ArcGeometry::SetCurvature(double curv)
{
	h = curv;
}

inline void ArcGeometry::SetTilt(double t)
{
	tilt = t;
}

inline double ArcGeometry::GetTilt() const
{
	return tilt;
}

#endif
