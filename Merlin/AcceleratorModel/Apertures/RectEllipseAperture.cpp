#include "AcceleratorModel/Apertures/RectEllipseAperture.h"
#include <iostream>
#include <cstdlib>
/*
* Checks if a particle is within a RectEllipseAperture type
* Returns true if the point (x,y,z) is within the aperture.
*/
bool RectEllipseAperture::PointInside (double x, double y, double z) const
{
	/* Particle is NOT inside the eliptical aperture component */
	//if( ((x*x*IEHH2) + (y*y*IEHV2)) > 1){return false;}
	if( (x*x + y*y*HV) > EHH2)
	{
		return false;
	}
	/* Particle is NOT inside the rectangular aperture component */
	else if(std::fabs(x) > RectHalfWidth || std::fabs(y) > RectHalfHeight)
	{
		return false;
	}
	/* Particle is inside both components, and is inside the aperture */
	else
	{
		return true;
	}
}

double RectEllipseAperture::GetRadiusAt (double phi, double z) const
{
	double a1 = RectHalfWidth;
	double a2 = RectHalfHeight;
	double a3 = EllipseHalfHorizontal;
	double a4 = EllipseHalfVertical;

	double t = phi;
	double rect_x = a1*((fabs(cos(t))*cos(t)) + (fabs(sin(t))*sin(t)));
	double rect_y = a2*((fabs(cos(t))*cos(t)) - (fabs(sin(t))*sin(t)));

	//Get the angle for the rectangle coordinate
	double EllipseAngle = atan2(rect_y,rect_x);
	double rr = a3*a4 / sqrt(pow(a4*cos(EllipseAngle),2) + pow(a3*sin(EllipseAngle),2));
	double ellipse_x = rr*cos(EllipseAngle);
	double ellipse_y = rr*sin(EllipseAngle);

	double rRectangle = sqrt((rect_x*rect_x) + (rect_y*rect_y));
	double rEllipse = sqrt((ellipse_x*ellipse_x) + (ellipse_y*ellipse_y));

	if(rRectangle < rEllipse)
	{
		return rRectangle;
	}
	else
	{
		return rEllipse;
	}
}

std::string RectEllipseAperture::GetApertureType() const
{
	return "RECTELLIPSE";
}

void RectEllipseAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(" << RectHalfWidth << ", "<< RectHalfHeight << ", "<<
	    EllipseHalfHorizontal << ", "<<EllipseHalfVertical << ")";
}

