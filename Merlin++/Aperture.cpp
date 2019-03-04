/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

#include "Aperture.h"
#include "CollimatorAperture.h"
#include "NumericalConstants.h"

using namespace std;

Aperture::Aperture()
{
	//for collimator inheritence
}

Aperture::~Aperture()
{
	//for collimator inheritence
}

Aperture::Aperture(string type, double s, double aper1) :
	apType(type), latticelocation(s)
{

}

Aperture::Aperture(string type, double s, double aper1, double aper2) :
	apType(type), latticelocation(s)
{

}

Aperture::Aperture(string type, double s, double aper1, double aper2, double aper3, double aper4) :
	apType(type), latticelocation(s)
{

}

void Aperture::printout(std::ostream& out) const
{
	out << apType;
}

CircularAperture::CircularAperture(double radius) :
	radius(radius)
{
	radius_sq = radius * radius;
}

CircularAperture::CircularAperture(string type, double s, double radius) :
	Aperture(type, s, radius), radius(radius)
{
	radius_sq = radius * radius;
}

string CircularAperture::GetType()
{
	return "CIRCLE";

}

Aperture* CircularAperture::GetInstance(DataTableRow dt)
{
	//THIS SHOULD BE APER_3 ACCORDING TO MADX USER GUIDE, BUT
	//CIRCLE APERTURE TWISS OUTPUT SEEMS TO PUT IT IN APER_1
	if(dt.Get_d("APER_3") == 0)
		return new CircularAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_1"));
	return new CircularAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_3"));
}

inline bool CircularAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < radius)
		return true;
	if(x * x + y * y >= radius_sq)
		return false;
	else
		return true;
}

void CircularAperture::printout(std::ostream& out) const
{
	out << apType << "(" << radius << ")";
}

RectangularAperture::RectangularAperture(double rectHalfX, double rectHalfY) :
	rectHalfX(rectHalfX), rectHalfY(rectHalfY)
{
	minDim = fmin(rectHalfX, rectHalfY);
	maxDim = fmax(rectHalfX, rectHalfY);
}

RectangularAperture::RectangularAperture(string type, double s, double rectHalfX, double rectHalfY) :
	Aperture(type, s, rectHalfX, rectHalfY), rectHalfX(rectHalfX), rectHalfY(rectHalfY)
{
	minDim = fmin(rectHalfX, rectHalfY);
	maxDim = fmax(rectHalfX, rectHalfY);
}

string RectangularAperture::GetType()
{
	return "RECTANGLE";

}

Aperture* RectangularAperture::GetInstance(DataTableRow dt)
{
	return new RectangularAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_1"), dt.Get_d("APER_2"));
}

inline bool RectangularAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minDim)
		return true;
	if(ax >= rectHalfX || ay >= rectHalfY)
		return false;
	else
		return true;
}

void RectangularAperture::printout(std::ostream& out) const
{
	out << apType << "(" << rectHalfX << ", " << rectHalfY << ")";
}

EllipticalAperture::EllipticalAperture(string type, double s, double ellipHalfX, double ellipHalfY) :
	Aperture(type, s, ellipHalfX, ellipHalfY), ellipHalfX(ellipHalfX), ellipHalfY(ellipHalfY)
{
	minDim = fmin(ellipHalfX, ellipHalfY);
	maxDim = fmax(ellipHalfX, ellipHalfY);
	ellipHalfX_sq = ellipHalfX * ellipHalfX;
	ellipHalfY_sq = ellipHalfY * ellipHalfY;
	ellipHalfX_sq_div_ellipHalfY_sq = ellipHalfX_sq / ellipHalfY_sq;
}

string EllipticalAperture::GetType()
{
	return "ELLIPSE";

}

Aperture* EllipticalAperture::GetInstance(DataTableRow dt)
{
	//SAME MAD-X TWISS OUT CONFUSION
	if(dt.Get_d("APER_3") == 0 && dt.Get_d("APER_4") == 0)
		return new EllipticalAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_1"), dt.Get_d("APER_2"));
	return new EllipticalAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_3"), dt.Get_d("APER_4"));
}

inline bool EllipticalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minDim)
		return true;
	if((x * x + y * y * ellipHalfX_sq_div_ellipHalfY_sq) >= ellipHalfX_sq)
		return false;
	else
		return true;
}

void EllipticalAperture::printout(std::ostream& out) const
{
	out << apType << "(" << ellipHalfX << ", " << ellipHalfY << ")";
}

RectEllipseAperture::RectEllipseAperture()
{
	//for collimator inheritence
}

RectEllipseAperture::RectEllipseAperture(string type, double s, double rectHalfX, double rectHalfY, double ellipHalfX,
	double
	ellipHalfY) :
	Aperture(type, s, rectHalfX, rectHalfY, ellipHalfX, ellipHalfY), rectHalfX(rectHalfX), rectHalfY(rectHalfY),
	ellipHalfX(ellipHalfX), ellipHalfY(ellipHalfY)
{
	minDim = min({rectHalfX, rectHalfY, ellipHalfX, ellipHalfY});
	maxDim = max({rectHalfX, rectHalfY, ellipHalfX, ellipHalfY});
	minRectDim = fmin(rectHalfX, rectHalfY);
	maxRectDim = fmax(rectHalfX, rectHalfY);
	minEllipDim = fmin(ellipHalfX, ellipHalfY);
	maxEllipDim = fmax(ellipHalfX, ellipHalfY);
	ellipHalfX_sq = ellipHalfX * ellipHalfX;
	ellipHalfY_sq = ellipHalfY * ellipHalfY;
	ellipHalfX_sq_div_ellipHalfY_sq = ellipHalfX_sq / ellipHalfY_sq;
}

string RectEllipseAperture::GetType()
{
	return "RECTELLIPSE";

}

Aperture* RectEllipseAperture::GetInstance(DataTableRow dt)
{
	return new RectEllipseAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_1"), dt.Get_d("APER_2"),
			   dt.Get_d("APER_3"), dt.Get_d("APER_4"));

}

inline bool RectEllipseAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minDim)
		return true;
	if((x * x + y * y * ellipHalfX_sq_div_ellipHalfY_sq) >= ellipHalfX_sq)
		return false;
	if(ax >= rectHalfX || ay >= rectHalfY)
		return false;
	else
		return true;
}

void RectEllipseAperture::printout(std::ostream& out) const
{
	out << apType << "(" << rectHalfX << ", " << rectHalfY << ", " << ellipHalfX << ", " << ellipHalfY << ")";
}

OctagonalAperture::OctagonalAperture(string type, double s, double rectHalfX, double rectHalfY, double angle1, double
	angle2) :
	Aperture(type, s, rectHalfX, rectHalfY, angle1, angle2), rectHalfX(rectHalfX), rectHalfY(rectHalfY), angle1(angle1),
	angle2(angle2)
{
	minDim = fmin(rectHalfX, rectHalfY);
	maxDim = fmax(rectHalfX, rectHalfY);
	const1 = (rectHalfY * tan((pi / 2) - angle2) - rectHalfX);
	const2 = rectHalfX * tan(angle1);
	const3 = rectHalfY - const2;
}

string OctagonalAperture::GetType()
{
	return "OCTAGON";

}

Aperture* OctagonalAperture::GetInstance(DataTableRow dt)
{
	return new OctagonalAperture(dt.Get_s("APERTYPE"), dt.Get_d("S"), dt.Get_d("APER_1"), dt.Get_d("APER_2"), dt.Get_d(
				   "APER_3"), dt.Get_d("APER_4"));

}

inline bool OctagonalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	//based on algorithm in trrun.f90 in MAD-X
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minDim)
		return true;
	if(ax >= rectHalfX || ay >= rectHalfY)
		return false;
	if(const1 * (ay - const2) - const3 * (ax - rectHalfX) <= 0)
		return false;
	return true;
}

void OctagonalAperture::printout(std::ostream& out) const
{
	out << apType << "(" << rectHalfX << ", " << rectHalfY << ", " << angle1 << ", " << angle2 << ")";
}

Aperture* ApertureFactory::GetInstance(DataTableRow dt)
{
	map<string, getAperture>::iterator itr = ApertureTypes.find(dt.Get_s("APERTYPE"));
	if(itr != ApertureTypes.end())
	{
		return (*itr->second)(dt);
	}
	return NULL;
}

ApertureFactoryInitializer::ApertureFactoryInitializer()
{
	ApertureFactory::ApertureTypes["CIRCLE"] = &CircularAperture::GetInstance;
	ApertureFactory::ApertureTypes["RECTANGLE"] = &RectangularAperture::GetInstance;
	ApertureFactory::ApertureTypes["ELLIPSE"] = &EllipticalAperture::GetInstance;
	ApertureFactory::ApertureTypes["RECTELLIPSE"] = &RectEllipseAperture::GetInstance;
	ApertureFactory::ApertureTypes["OCTAGON"] = &OctagonalAperture::GetInstance;
}

map<string, getAperture> ApertureFactory::ApertureTypes;
ApertureFactoryInitializer ApertureFactoryInitializer::init;
