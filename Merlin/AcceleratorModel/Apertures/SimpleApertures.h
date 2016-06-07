/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
//
/////////////////////////////////////////////////////////////////////////

#ifndef SimpleApertures_h
#define SimpleApertures_h 1

#include "merlin_config.h"
#include <cmath>
#include "AcceleratorModel/Aperture.h"
#include "NumericalUtils/NumericalConstants.h"

#include <string>
#include <iostream>

#ifdef MERLIN_PROFILE
#include "utility/MerlinProfile.h"
#endif

/**
* Represents an aperture with a rectangular cross-section.
* The aperture is assumed symmetric about the axis, and
* extruded along its geometry.
*/
class RectangularAperture : public Aperture
{
public:
	RectangularAperture (double width, double height);

	double GetFullWidth () const;
	double GetFullHeight () const;
	void SetFullWidth (double w);
	void SetFullHeight (double h);

	/**
	* Returns true if the point (x,y,z) is within the
	* aperture. The z coordinate is ignored.
	*/
	virtual bool PointInside (double x, double y, double z) const;

	/**
	* Returns the radius to the aperture at the angle phi. The
	* z coordinate is ignored.
	*/
	virtual double GetRadiusAt (double phi, double z) const;

	virtual std::string GetApertureType() const;
	virtual void printout(std::ostream& out) const;

private:

	double hw;
	double hh;
};

inline std::string RectangularAperture::GetApertureType() const
{
	return "RECTANGULAR";
}

/**
* Represents an aperture with a circular cross-section.
* The aperture is assumed to be extruded along its
* geometry.
*/
class CircularAperture : public Aperture
{
public:
	explicit CircularAperture (double r);

	double GetRadius () const;
	double GetDiameter () const;
	void SetRadius (double r);
	void SetDiameter (double d);

	/**
	* Returns true if the point (x,y,z) is within the
	* aperture. The z coordinate is ignored.
	*/
	virtual bool PointInside (double x, double y, double z) const;

	/**
	* Returns the radius.
	*/
	virtual double GetRadiusAt (double phi, double z) const;

	virtual std::string GetApertureType() const;
	virtual void printout(std::ostream& out) const;
private:

	double r2;
};

inline RectangularAperture::RectangularAperture (double width, double height)
	:hw(fabs(width)/2),hh(fabs(height)/2)
{}

inline double RectangularAperture::GetFullWidth () const
{
	return 2*hw;
}

inline double RectangularAperture::GetFullHeight () const
{
	return 2*hh;
}

inline void RectangularAperture::SetFullWidth (double w)
{
	hw=fabs(w)/2;
}

inline void RectangularAperture::SetFullHeight (double h)
{
	hh=fabs(h)/2;
}

inline bool RectangularAperture::PointInside (double x, double y, double z=0.0) const
{
	return fabs(x)<hw && fabs(y)<hh;
}

inline CircularAperture::CircularAperture (double r)
	: r2(r*r)
{}

inline double CircularAperture::GetRadius () const
{
	return sqrt(r2);
}

inline double CircularAperture::GetDiameter () const
{
	return 2*GetRadius();
}

inline void CircularAperture::SetRadius (double r)
{
	r2=r*r;
}

inline void CircularAperture::SetDiameter (double d)
{
	r2=d*d/4.0;
}

inline bool CircularAperture::PointInside (double x, double y, double z) const
{
	return x*x+y*y<r2;
}

inline std::string CircularAperture::GetApertureType() const
{
	return "CIRCULAR";
}

/**
* Represents an aperture with an elliptical cross-section.
* The aperture is assumed extruded along its geometry.
*/
class EllipticalAperture : public Aperture
{
public:
	EllipticalAperture (double width, double height);

	double GetHalfWidth () const;
	double GetHalfHeight () const;

	/**
	* Returns true if the point (x,y,z) is within the
	* aperture. The z coordinate is ignored.
	*/
	virtual bool PointInside (double x, double y, double z) const;

	/**
	* Returns the radius to the aperture at the angle phi. The
	* z coordinate is ignored.
	*/
	virtual double GetRadiusAt (double phi, double z) const;

	virtual std::string GetApertureType() const;
	virtual void printout(std::ostream& out) const;

private:

	double hw;
	double hh;
	double EHH2;
	double HV;
};

inline EllipticalAperture::EllipticalAperture (double ehh, double ehv)
	: hw(ehh), hh(ehv), EHH2(ehh*ehh), HV((ehh*ehh)/(ehv*ehv))
{}

inline double EllipticalAperture::GetHalfWidth () const
{
	return hw;
}

inline double EllipticalAperture::GetHalfHeight () const
{
	return hh;
}

inline bool EllipticalAperture::PointInside (double x, double y, double z) const
{
	return (x*x + y*y*HV) < EHH2;
}

inline std::string EllipticalAperture::GetApertureType() const
{
	return "ELLIPSE";
}

/**
* Represents an aperture with an octagon cross-section.
* The aperture is assumed extruded along its geometry.
*/
class OctagonalAperture : public Aperture
{
public:
	OctagonalAperture (double width, double height, double ang1, double ang2);

	//Functions to extract the aperture parameters
	double GetHalfWidth() const;
	double GetHalfHeight() const;
	double GetAngle1() const;
	double GetAngle2() const;

	/**
	* Returns true if the point (x,y,z) is within the
	* aperture. The z coordinate is ignored.
	*/
	virtual bool PointInside (double x, double y, double z) const;

	/**
	* Returns the radius to the aperture at the angle phi. The
	* z coordinate is ignored.
	*/
	virtual double GetRadiusAt (double phi, double z) const;

	virtual std::string GetApertureType() const;
	virtual void printout(std::ostream& out) const;

private:

	//Aperture parameters
	double hw;
	double hh;
	double angle1;
	double angle2;

	//This tangents are frequently used and are the same each time so will be pre-computed
	double tana1;
	double tana2;

	double c1;
	double c2;
	double c3;
};

inline OctagonalAperture::OctagonalAperture (double h, double v, double a1, double a2)
	: hw(h), hh(v), angle1(a1), angle2(a2)
{
	//Compute the tangents.
	tana1 = tan(angle1);
	tana2 = tan(pi/2 - angle2);

	//Make some constants.
	//(hh*tana2 - hw)
	c1 = (hh*tana2 - hw);

	//hw*tana1
	c2 = hw*tana1;

	//hh - hw*tana1
	c3 = hh - c2;
}

inline bool OctagonalAperture::PointInside (double x, double y, double z) const
{
	//This is just taken from trrun.f90 in MAD-X. - credit to: 2015-Feb-20  18:42:26  ghislain: added octagon shape

	/*
	!*** case of octagon: test outer rectangle (ap1,ap2) then test cut corner.
	lost =  x .gt. ap1 .or. y .gt. ap2 .or. &
	     (ap2*tan(pi/2 - ap4) - ap1)*(y - ap1*tan(ap3)) - (ap2 - ap1*tan(ap3))*(x - ap1) .lt. zero
	*/

	double fabsx = fabs(x);
	double fabsy = fabs(y);

	x=fabsx;
	y=fabsy;
	//First check the rectangle
	if(x >= hw || y >= hh)
	{
		return false;
	}

	if(c1*(y - c2) - c3*(x - hw) <= 0 )
	{
		return false;
	}

	//Particle survives both checks
	return true;
}

inline double OctagonalAperture::GetHalfWidth() const
{
	return hw;
}

inline double OctagonalAperture::GetHalfHeight() const
{
	return hh;
}

inline double OctagonalAperture::GetAngle1() const
{
	return angle1;
}

inline double OctagonalAperture::GetAngle2() const
{
	return angle2;
}

inline std::string OctagonalAperture::GetApertureType() const
{
	return "OCTAGON";
}

#endif

