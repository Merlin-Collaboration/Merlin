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

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "NumericalUtils/NumericalConstants.h"
#include "NumericalUtils/utils.h"

#include <iostream>
#include <cstdlib>

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

double CircularAperture::GetRadiusAt (double phi, double z) const
{
	return GetRadius();
}

void CircularAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << GetRadius () <<")";
}

