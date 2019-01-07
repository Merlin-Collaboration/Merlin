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

using namespace std;

ApertureAbstract::ApertureAbstract()
{

}

ApertureAbstract::~ApertureAbstract()
{

}

Aperture::Aperture()
{

}

Aperture::Aperture(string type, double s, double aper1, double aper2, double aper3, double aper4) :
	apType(type), s_longitudinal(s), rectHalfWidth(aper1), rectHalfHeight(aper2), ellipHalfWidth(aper3),
	ellipHalfHeight(aper4), minRectDim(0), minEllipDim(0), minDim(0), maxRectDim(0), ellipHalfHeight2(0),
	ellipHalfWidth2(0), ellipHalfWidth2overEllipHalfHeight2(0)
{
	CalcApertureParams(rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight);
}

CircularAperture::CircularAperture(double aper3)
{
	setEllipHalfWidth(aper3);
}

CircularAperture::CircularAperture(string type, double s, double aper1, double aper2, double aper3, double aper4) :
	Aperture(type, s, aper1, aper2, aper3, aper4)
{

}

string CircularAperture::getType()
{
	return "CIRCLE";

}

Aperture* CircularAperture::getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4)
{
	return new CircularAperture(type, s, aper1, aper2, aper3, aper4);

}

inline bool CircularAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minEllipDim)
		return true;
	if(x * x + y * y >= ellipHalfWidth2)
		return false;
	else
		return true;
}

RectangularAperture::RectangularAperture(double aper1, double aper2)
{
	setRectHalfWidth(aper1);
	setRectHalfHeight(aper2);
}

RectangularAperture::RectangularAperture(string type, double s, double aper1, double aper2, double aper3, double
	aper4) :
	Aperture(type, s, aper1, aper2, aper3, aper4)
{

}

string RectangularAperture::getType()
{
	return "RECTANGLE";

}

Aperture* RectangularAperture::getInstance(string type, double s, double aper1, double aper2, double aper3, double
	aper4)
{
	return new RectangularAperture(type, s, aper1, aper2, aper3, aper4);

}

inline bool RectangularAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minRectDim)
		return true;
	if(ax >= rectHalfWidth || ay >= rectHalfHeight)
		return false;
	else
		return true;
}

EllipticalAperture::EllipticalAperture(double aper3, double aper4)
{
	setEllipHalfWidth(aper3);
	setEllipHalfHeight(aper4);
}

EllipticalAperture::EllipticalAperture(string type, double s, double aper1, double aper2, double aper3, double aper4) :
	Aperture(type, s, aper1, aper2, aper3, aper4)
{

}

string EllipticalAperture::getType()
{
	return "ELLIPSE";

}

Aperture* EllipticalAperture::getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4)
{
	return new EllipticalAperture(type, s, aper1, aper2, aper3, aper4);

}

inline bool EllipticalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minEllipDim)
		return true;
	if((x * x + y * y * ellipHalfWidth2overEllipHalfHeight2) >= ellipHalfWidth2)
		return false;
	else
		return true;
}

RectEllipseAperture::RectEllipseAperture()
{

}

RectEllipseAperture::RectEllipseAperture(double aper1, double aper2, double aper3, double aper4)
{
	setRectHalfWidth(aper1);
	setRectHalfHeight(aper2);
	setEllipHalfWidth(aper3);
	setEllipHalfHeight(aper4);
}

RectEllipseAperture::RectEllipseAperture(string type, double s, double aper1, double aper2, double aper3, double
	aper4) :
	Aperture(type, s, aper1, aper2, aper3, aper4)
{

}

string RectEllipseAperture::getType()
{
	return "RECTELLIPSE";

}

Aperture* RectEllipseAperture::getInstance(string type, double s, double aper1, double aper2, double aper3, double
	aper4)
{
	return new RectEllipseAperture(type, s, aper1, aper2, aper3, aper4);

}

inline bool RectEllipseAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minDim)
		return true;
	if((x * x + y * y * ellipHalfWidth2overEllipHalfHeight2) >= ellipHalfWidth2)
		return false;
	if(ax >= rectHalfWidth || ay >= rectHalfHeight)
		return false;
	else
		return true;
}

OctagonalAperture::OctagonalAperture()
{

}

OctagonalAperture::OctagonalAperture(string type, double s, double aper1, double aper2, double aper3, double aper4) :
	Aperture(type, s, aper1, aper2, aper3, aper4), angle1(aper3), angle2(aper4), const1(
		rectHalfHeight * tan(angle1) - rectHalfWidth), const2(rectHalfWidth * tan(angle1)), const3(
		rectHalfHeight - const2)
{

}

string OctagonalAperture::getType()
{
	return "OCTAGON";

}

Aperture* OctagonalAperture::getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4)
{
	return new OctagonalAperture(type, s, aper1, aper2, aper3, aper4);

}

inline bool OctagonalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double ax = fabs(x);
	double ay = fabs(y);
	if(ax + ay < minRectDim)
		return true;
	if(ax >= rectHalfWidth || ay >= rectHalfHeight)
		return false;
	if(const1 * (y - const2) - const3 * (x - rectHalfWidth) <= 0)
		return false;
	else
		return true;
}

Aperture* ApertureFactory::getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4)
{
	map<string, getAperture>::iterator itr = apertureTypes.find(type);
	if(itr != apertureTypes.end())
	{
		return (*itr->second)(type, s, aper1, aper2, aper3, aper4);
	}
	return NULL;
}

ApertureFactoryInitializer::ApertureFactoryInitializer()
{
	ApertureFactory::apertureTypes["CIRCLE"] = &CircularAperture::getInstance;
	ApertureFactory::apertureTypes["RECTANGLE"] = &RectangularAperture::getInstance;
	ApertureFactory::apertureTypes["ELLIPSE"] = &EllipticalAperture::getInstance;
	ApertureFactory::apertureTypes["RECTELLIPSE"] = &RectEllipseAperture::getInstance;
	ApertureFactory::apertureTypes["OCTAGON"] = &OctagonalAperture::getInstance;
}

map<string, getAperture> ApertureFactory::apertureTypes;
ApertureFactoryInitializer ApertureFactoryInitializer::init;
