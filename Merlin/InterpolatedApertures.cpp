#include "InterpolatedApertures.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

/**
* InterpolatedRectEllipseAperture
*/

//Returns true if the point (x,y,z) is within the aperture.
bool InterpolatedRectEllipseAperture::PointInside (double x, double y, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack = InterpolatedAperture::ApertureList[n-1];
			apFront = InterpolatedAperture::ApertureList[n];
			break;
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g1 = (apFront.ap1 - apBack.ap1) / delta_s;
	double g2 = (apFront.ap2 - apBack.ap2) / delta_s;
	double g3 = (apFront.ap3 - apBack.ap3) / delta_s;
	double g4 = (apFront.ap4 - apBack.ap4) / delta_s;

	//c = y - mx
	double c1 = apFront.ap1 - (g1 * apFront.s);
	double c2 = apFront.ap2 - (g2 * apFront.s);
	double c3 = apFront.ap3 - (g3 * apFront.s);
	double c4 = apFront.ap4 - (g4 * apFront.s);

	double rect_half_width = (g1 * z) + c1;
	double rect_half_height = (g2 * z) + c2;
	double ellipse_half_horizontal = (g3 * z) + c3;
	double ellipse_half_vertical = (g4 * z) + c4;

	//aper_1 = half width rectangle
	//aper_2 = half heigth rectangle
	//aper_3 = half horizontal axis ellipse (or radius if circle)
	//aper_4 = half vertical axis ellipse

	if(((x*x)/(ellipse_half_horizontal*ellipse_half_horizontal)) + ((y*y)/(ellipse_half_vertical*ellipse_half_vertical)) > 1)
	{
		return 0;
	}
	else if(fabs(x) > rect_half_width || fabs(y) > rect_half_height)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double InterpolatedRectEllipseAperture::GetRadiusAt (double phi, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack = InterpolatedAperture::ApertureList[n-1];
			apFront = InterpolatedAperture::ApertureList[n];
			break;
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g1 = (apFront.ap1 - apBack.ap1) / delta_s;
	double g2 = (apFront.ap2 - apBack.ap2) / delta_s;
	double g3 = (apFront.ap3 - apBack.ap3) / delta_s;
	double g4 = (apFront.ap4 - apBack.ap4) / delta_s;

	//c = y - mx
	double c1 = apFront.ap1 - (g1 * apFront.s);
	double c2 = apFront.ap2 - (g2 * apFront.s);
	double c3 = apFront.ap3 - (g3 * apFront.s);
	double c4 = apFront.ap4 - (g4 * apFront.s);

	double rect_half_width = (g1 * z) + c1;
	double rect_half_height = (g2 * z) + c2;
	double ellipse_half_horizontal = (g3 * z) + c3;
	double ellipse_half_vertical = (g4 * z) + c4;

	double a1 = rect_half_width;
	double a2 = rect_half_height;
	double a3 = ellipse_half_horizontal;
	double a4 = ellipse_half_vertical;

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


void InterpolatedRectEllipseAperture::EnablePrint()
{
	Print = true;
}

void InterpolatedRectEllipseAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(";
	for(size_t n=0; n < ApertureList.size(); n++)
	{
		out << ApertureList[n].s << " [";
		out << ApertureList[n].ap1<< ", " << ApertureList[n].ap2<< ", "<< ApertureList[n].ap3<< ", "<< ApertureList[n].ap4;
		out << "]";
		if (n < ApertureList.size()-1)
		{
			out << ", ";
		}
	}
	out << ")";
}

/**
* InterpolatedCircularAperture
*/

bool InterpolatedCircularAperture::PointInside (double x, double y, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack.s = InterpolatedAperture::ApertureList[n-1].s;
			apBack.ap1 = InterpolatedAperture::ApertureList[n-1].ap1;

			apFront.s = InterpolatedAperture::ApertureList[n].s;
			apFront.ap1 = InterpolatedAperture::ApertureList[n].ap1;

			if(apFront.ap1 == 0 || apBack.ap1 == 0)
			{
				std::cout << z << std::endl;
				std::cout << apFront.s << std::endl;
				std::cout << apBack.s << std::endl;
				std::cout << apFront.ap1 << std::endl;
				std::cout << apBack.ap1 << std::endl;
				std::cout << ApertureList[n].ap1 << std::endl;
				std::cout << ApertureList[n-1].ap1 << std::endl;
				abort();
			}

			break;
		}

		if(n == (ApertureList.size()-1))
		{
			std::cout << "No aperture found for InterpolatedCircularAperture at z = " << z << std::endl;
			for(size_t m=0; m < ApertureList.size(); m++)
			{
				std::cout << "Entry " << m << " - z = " << ApertureList[m].s << "\t" << ApertureList[m].ap1 << std::endl;
			}
			abort();
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g = (apFront.ap1 - apBack.ap1) / delta_s;

	//c = y - mx
	double c = apFront.ap1 - (g * apFront.s);

	//aper_3 = half horizontal axis ellipse (or radius if circle)
	double r2 = pow((g * z) + c,2);

	if((x*x+y*y<r2) == 0)
	{
		return false;
	}

	return x*x+y*y<r2;
}

double InterpolatedCircularAperture::GetRadiusAt (double phi, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack = InterpolatedAperture::ApertureList[n-1];
			apFront = InterpolatedAperture::ApertureList[n];
			break;
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g = (apFront.ap1 - apBack.ap1) / delta_s;

	//c = y - mx
	double c = apFront.ap1 - (g * apFront.s);

	//aper_3 = half horizontal axis ellipse (or radius if circle)
	double r2 = pow((g * z) + c,2);
	return sqrt(r2);
}

double InterpolatedCircularAperture::GetRadius () const
{
	return GetRadiusAt(0.0, 0.0);
}

void InterpolatedCircularAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(";
	for(size_t n=0; n < ApertureList.size(); n++)
	{
		out << ApertureList[n].s << " [";
		out << ApertureList[n].ap1;
		out << "]";
		if (n < ApertureList.size()-1)
		{
			out << ", ";
		}
	}
	out << ")";
}

/**
* InterpolatedEllipticalAperture
*/

bool InterpolatedEllipticalAperture::PointInside (double x, double y, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack.s = InterpolatedAperture::ApertureList[n-1].s;
			apBack.ap1 = InterpolatedAperture::ApertureList[n-1].ap1;

			apFront.s = InterpolatedAperture::ApertureList[n].s;
			apFront.ap1 = InterpolatedAperture::ApertureList[n].ap1;

			if(apFront.ap1 == 0 || apBack.ap1 == 0)
			{
				std::cout << z << std::endl;
				std::cout << apFront.s << std::endl;
				std::cout << apBack.s << std::endl;
				std::cout << apFront.ap1 << std::endl;
				std::cout << apBack.ap1 << std::endl;
				std::cout << ApertureList[n].ap1 << std::endl;
				std::cout << ApertureList[n-1].ap1 << std::endl;
				abort();
			}

			break;
		}

		if(n == (ApertureList.size()-1))
		{
			std::cout << "No aperture found for InterpolatedEllipticalAperture at z = " << z << std::endl;
			for(size_t m=0; m < ApertureList.size(); m++)
			{
				std::cout << "Entry " << m << " - z = " << ApertureList[m].s << "\t" << ApertureList[m].ap1 << std::endl;
			}
			abort();
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g1 = (apFront.ap1 - apBack.ap1) / delta_s;
	double g2 = (apFront.ap2 - apBack.ap2) / delta_s;

	//c = y - mx
	double c1 = apFront.ap1 - (g1 * apFront.s);
	double c2 = apFront.ap1 - (g2 * apFront.s);

	double ellipse_half_horizontal = (g1 * z) + c1;
	double ellipse_half_vertical = (g2 * z) + c2;

	if(((x*x)/(ellipse_half_horizontal*ellipse_half_horizontal)) + ((y*y)/(ellipse_half_vertical*ellipse_half_vertical)) > 1)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double InterpolatedEllipticalAperture::GetRadiusAt (double phi, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack.s = InterpolatedAperture::ApertureList[n-1].s;
			apBack.ap1 = InterpolatedAperture::ApertureList[n-1].ap1;

			apFront.s = InterpolatedAperture::ApertureList[n].s;
			apFront.ap1 = InterpolatedAperture::ApertureList[n].ap1;

			if(apFront.ap1 == 0 || apBack.ap1 == 0)
			{
				std::cout << z << std::endl;
				std::cout << apFront.s << std::endl;
				std::cout << apBack.s << std::endl;
				std::cout << apFront.ap1 << std::endl;
				std::cout << apBack.ap1 << std::endl;
				std::cout << ApertureList[n].ap1 << std::endl;
				std::cout << ApertureList[n-1].ap1 << std::endl;
				abort();
			}

			break;
		}

		if(n == (ApertureList.size()-1))
		{
			std::cout << "No aperture found for InterpolatedEllipticalAperture at z = " << z << std::endl;
			for(size_t m=0; m < ApertureList.size(); m++)
			{
				std::cout << "Entry " << m << " - z = " << ApertureList[m].s << "\t" << ApertureList[m].ap1 << std::endl;
			}
			abort();
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g1 = (apFront.ap1 - apBack.ap1) / delta_s;
	double g2 = (apFront.ap2 - apBack.ap2) / delta_s;

	//c = y - mx
	double c1 = apFront.ap1 - (g1 * apFront.s);
	double c2 = apFront.ap1 - (g2 * apFront.s);

	double ellipse_half_horizontal = (g1 * z) + c1;
	double ellipse_half_vertical = (g2 * z) + c2;

	double rr = ellipse_half_horizontal*ellipse_half_vertical / sqrt(pow(ellipse_half_vertical*cos(phi),2) + pow(ellipse_half_horizontal*sin(phi),2));
	double ellipse_x = rr*cos(phi);
	double ellipse_y = rr*sin(phi);
	return sqrt((ellipse_x*ellipse_x) + (ellipse_y*ellipse_y));
}

void InterpolatedEllipticalAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(";
	for(size_t n=0; n < ApertureList.size(); n++)
	{
		out << ApertureList[n].s << " [";
		out << ApertureList[n].ap1<< ", " << ApertureList[n].ap2;
		out << "]";
		if (n < ApertureList.size()-1)
		{
			out << ", ";
		}
	}
	out << ")";
}

/**
* InterpolatedOctagonalAperture
*/

//Returns true if the point (x,y,z) is within the aperture.
bool InterpolatedOctagonalAperture::PointInside (double x, double y, double z) const
{
	ap apFront;
	ap apBack;

	apFront.s = 0;
	apFront.ap1 = 0;
	apFront.ap2 = 0;
	apFront.ap3 = 0;
	apFront.ap4 = 0;

	apBack.s = 0;
	apBack.ap1 = 0;
	apBack.ap2 = 0;
	apBack.ap3 = 0;
	apBack.ap4 = 0;

	if(z < 0)
	{
		z = 0;
	}

	//The z coordinate is used to find the other aperture components
	for(size_t n=1; n < ApertureList.size(); n++)
	{
		if(ApertureList[n].s >= z)
		{
			apBack = InterpolatedAperture::ApertureList[n-1];
			apFront = InterpolatedAperture::ApertureList[n];
			break;
		}
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently.

	//y = mx + c
	//m = (y1 - y0) / (x1 - x0)
	double delta_s = apFront.s - apBack.s;
	double g1 = (apFront.ap1 - apBack.ap1) / delta_s;
	double g2 = (apFront.ap2 - apBack.ap2) / delta_s;
	double g3 = (apFront.ap3 - apBack.ap3) / delta_s;
	double g4 = (apFront.ap4 - apBack.ap4) / delta_s;

	//c = y - mx
	double c1 = apFront.ap1 - (g1 * apFront.s);
	double c2 = apFront.ap2 - (g2 * apFront.s);
	double c3 = apFront.ap3 - (g3 * apFront.s);
	double c4 = apFront.ap4 - (g4 * apFront.s);

	// hw(h), hh(v), angle1(a1), angle2(a2)
	double hw = (g1 * z) + c1;
	double hh = (g2 * z) + c2;
	double angle1 = (g3 * z) + c3;
	double angle2 = (g4 * z) + c4;

	//Compute the tangents.
	double tana1 = tan(angle1);
	double tana2 = tan(pi/2 - angle2);

	//Make some constants.
	//(hh*tana2 - hw)
	double cc1 = (hh*tana2 - hw);

	//hw*tana1
	double cc2 = hw*tana1;

	//hh - hw*tana1
	double cc3 = hh - cc2;

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

	if(cc1*(y - cc2) - cc3*(x - hw) <= 0 )
	{
		return false;
	}

	//Particle survives both checks
	return true;
}

double InterpolatedOctagonalAperture::GetRadiusAt (double phi, double z) const
{
	std::cerr << "InterpolatedOctagonalAperture::GetRadiusAt() - Not yet implemented." << std::endl;
	exit(EXIT_FAILURE);
	return 0;
}


void InterpolatedOctagonalAperture::EnablePrint()
{
	Print = true;
}

void InterpolatedOctagonalAperture::printout(std::ostream& out) const
{
	out << GetApertureType() << "(";
	for(size_t n=0; n < ApertureList.size(); n++)
	{
		out << ApertureList[n].s << " [";
		out << ApertureList[n].ap1<< ", " << ApertureList[n].ap2<< ", "<< ApertureList[n].ap3<< ", "<< ApertureList[n].ap4;
		out << "]";
		if (n < ApertureList.size()-1)
		{
			out << ", ";
		}
	}
	out << ")";
}

