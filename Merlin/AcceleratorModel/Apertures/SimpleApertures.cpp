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
        return 0;

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

/*
//	Returns true if the point (x,y,z) is within the aperture.
inline bool RectEllipseAperture::PointInside (double x, double y, double z) const
{
	#ifdef MERLIN_PROFILE
	MerlinProfile::StartProcessTimer("APERTURE");
	#endif
	//cout << "RECTELLIPSE" << endl;
	if(((x*x)/(ellipse_half_horizontal*ellipse_half_horizontal)) + ((y*y)/(ellipse_half_vertical*ellipse_half_vertical)) > 1)
	{
		//Particle is NOT inside the eliptical aperture component
//		cout << "Not in ellipse" << endl;
		#ifdef MERLIN_PROFILE
		MerlinProfile::EndProcessTimer("APERTURE");
		#endif
		return 0;
	}
	else if(fabs(x) > rect_half_width || fabs(y) > rect_half_height)
	{
		//Particle is NOT inside the rectangular aperture component
//		cout << "Not in rectangle:\tx: " << x << "\t" << rect_half_width << "\ty: " << y << "\t" << rect_half_height << endl;
		#ifdef MERLIN_PROFILE
		MerlinProfile::EndProcessTimer("APERTURE");
		#endif
		return 0;
	}
	else
	{
		//Particle is inside both components, and is inside the aperture
//		cout << "RECTELLIPSE APERTURE CHECK ok" << endl;
		#ifdef MERLIN_PROFILE
		MerlinProfile::EndProcessTimer("APERTURE");
		#endif
		return 1;
	}
}

double RectEllipseAperture::GetRadiusAt (double phi, double z) const
{

	std::cerr << "Not yet implemented" << std::endl;
	exit(EXIT_FAILURE);
}
*/
