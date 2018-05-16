/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "SimpleApertures.h"
#include "NumericalConstants.h"
#include "utils.h"

#include <iostream>
#include <cstdlib>

/**
* Rectangular Aperture Functions
*/

double RectangularAperture::GetRadiusAt (double phi, double z) const
{
	if(fequal(hh,0.0) || fequal(hw,0.0))
	{
		return 0;
	}

	const double phi0=atan(hh/hw);
	const double piOverTwo = pi/2.0;

	phi=fmod(phi,piOverTwo)*piOverTwo;

	return phi<phi0 ? hw/cos(phi) : hh/sin(phi);
}

void RectangularAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << hw << ", " << hh <<")";
}

/**
* Circular Aperture Functions
*/

double CircularAperture::GetRadiusAt (double phi, double z) const
{
	return GetRadius();
}

void CircularAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << GetRadius () <<")";
}

/**
* Elliptical Aperture Functions
*/

double EllipticalAperture::GetRadiusAt (double phi, double z) const
{
	double rr = hw*hh / sqrt(pow(hh*cos(phi),2) + pow(hw*sin(phi),2));
	double ellipse_x = rr*cos(phi);
	double ellipse_y = rr*sin(phi);
	return sqrt((ellipse_x*ellipse_x) + (ellipse_y*ellipse_y));
}

void EllipticalAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << GetHalfWidth() << ", " << GetHalfHeight() << ")";
}

/**
* Octagonal Aperture Functions
*/

double OctagonalAperture::GetRadiusAt (double phi, double z) const
{
	std::cerr << "OctagonalAperture::GetRadiusAt() - not yet implemented" << std::endl;
	exit(EXIT_FAILURE);
}

void OctagonalAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << GetHalfWidth() << ", " << GetHalfHeight() << ", " << GetAngle1() << ", " << GetAngle2() << ")";
}

