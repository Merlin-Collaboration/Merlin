/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "NumericalConstants.h"
#include "InterpolatedApertures.h"

using namespace std;

string InterpolatedAperture::GetType()
{
	return "RECTELLIPSE-interpolated";
}

InterpolatedAperture::InterpolatedAperture(DataTable dt)
{
	AperturesToInterpolate = dt;
	ConvertToStruct(AperturesToInterpolate);
}

InterpolatedAperture::~InterpolatedAperture()
{

}

Aperture* InterpolatedAperture::GetInstance(DataTable dt)
{
	return new InterpolatedAperture(dt);
}

bool InterpolatedAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	//if you're here you're using a yet unsupported interpolation
	cout << "Unsupported Aperture Interpolation... exiting." << endl;
	exit(EXIT_FAILURE);
	return true;
}

void InterpolatedAperture::ConvertToStruct(DataTable dt)
{
	for(size_t row = 0; row < dt.Length(); ++row)
	{
		ApEntry.s = dt.Get_d("S",row);
		ApEntry.ap1 = dt.Get_d("APER_1",row);
		ApEntry.ap2 = dt.Get_d("APER_2",row);
		ApEntry.ap3 = dt.Get_d("APER_3",row);
		ApEntry.ap4 = dt.Get_d("APER_4",row);
		ApList.push_back(ApEntry);
	}
}

InterpolatedRectEllipseAperture::InterpolatedRectEllipseAperture(DataTable dt):
				InterpolatedAperture(dt)
{

}

InterpolatedRectEllipseAperture::~InterpolatedRectEllipseAperture()
{

}

Aperture* InterpolatedRectEllipseAperture::GetInstance(DataTable dt)
{
	return new InterpolatedRectEllipseAperture(dt);
}

string InterpolatedRectEllipseAperture::GetType()
{
	return "ELLIPSE-interpolated";
}

bool InterpolatedRectEllipseAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	apStruct apFront;
	apStruct apBack;
	double ax = fabs(x);
	double ay = fabs(y);


	if(z < 0)
	{
		z = 0;
	}

	for(size_t n = 1; n < ApList.size(); ++n)
	{
		if(ApList[n].s >= z)
		{
			apBack = ApList[n - 1];
			apFront = ApList[n];
			break;
		}
	}

	double minDim = min({apFront.ap1, apFront.ap2, apFront.ap3, apFront.ap4, apFront.ap1, apBack.ap2, apBack.ap3, apBack.ap4});
	if(ax + ay < minDim)
	{
		return true;
	}

	//We now have the aperture before and after the particle.
	//It must now be calculated at the point where the particle is currently
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

	if(((x * x) / (ellipse_half_horizontal * ellipse_half_horizontal)) + ((y * y) / (ellipse_half_vertical
			* ellipse_half_vertical)) > 1)
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

InterpolatedCircularAperture::InterpolatedCircularAperture(DataTable dt):
				InterpolatedAperture(dt)
{

}

InterpolatedCircularAperture::~InterpolatedCircularAperture()
{

}

Aperture* InterpolatedCircularAperture::GetInstance(DataTable dt)
{
	return new InterpolatedCircularAperture(dt);
}

string InterpolatedCircularAperture::GetType()
{
	return "CIRCLE-interpolated";
}

bool InterpolatedCircularAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	apStruct apFront;
	apStruct apBack;
	double ax = fabs(x);
	double ay = fabs(y);

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
	for(size_t n = 1; n < ApList.size(); n++)
	{
		if(ApList[n].s >= z)
		{
			apBack.s = ApList[n - 1].s;
			apBack.ap1 = ApList[n - 1].ap1;

			apFront.s = ApList[n].s;
			apFront.ap1 = ApList[n].ap1;


			double minDim = min({apFront.ap1, apBack.ap1});
			if(ax + ay < minDim)
			{
				return true;
			}

			if(apFront.ap1 == 0 || apBack.ap1 == 0)
			{
				cout << z << endl;
				cout << apFront.s << endl;
				cout << apBack.s << endl;
				cout << apFront.ap1 << endl;
				cout << apBack.ap1 << endl;
				cout << ApList[n].ap1 << endl;
				cout << ApList[n - 1].ap1 << endl;
				abort();
			}

			break;
		}

		if(n == (ApList.size() - 1))
		{
			cout << "No aperture found for InterpolatedCircularAperture at z = " << z << endl;
			for(size_t m = 0; m < ApList.size(); m++)
			{
				cout << "Entry " << m << " - z = " << ApList[m].s << "\t" << ApList[m].ap1
						<< endl;
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
	double r2 = pow((g * z) + c, 2);

	if((x * x + y * y < r2) == 0)
	{
		return false;
	}

	return x * x + y * y < r2;
}

InterpolatedEllipticalAperture::InterpolatedEllipticalAperture(DataTable dt):
				InterpolatedAperture(dt)
{

}

InterpolatedEllipticalAperture::~InterpolatedEllipticalAperture()
{

}

Aperture* InterpolatedEllipticalAperture::GetInstance(DataTable dt)
{
	return new InterpolatedEllipticalAperture(dt);
}

string InterpolatedEllipticalAperture::GetType()
{
	return "CIRCLE-interpolated";
}

bool InterpolatedEllipticalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	apStruct apFront;
	apStruct apBack;
	double ax = fabs(x);
	double ay = fabs(y);

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
	for(size_t n = 1; n < ApList.size(); n++)
	{
		if(ApList[n].s >= z)
		{
			apBack.s = ApList[n - 1].s;
			apBack.ap1 = ApList[n - 1].ap1;

			apFront.s = ApList[n].s;
			apFront.ap1 = ApList[n].ap1;

			double minDim = min({apFront.ap1, apBack.ap1});
			if(ax + ay < minDim)
			{
				return true;
			}

			if(apFront.ap1 == 0 || apBack.ap1 == 0)
			{
				cout << z << endl;
				cout << apFront.s << endl;
				cout << apBack.s << endl;
				cout << apFront.ap1 << endl;
				cout << apBack.ap1 << endl;
				cout << ApList[n].ap1 << endl;
				cout << ApList[n - 1].ap1 << endl;
				abort();
			}

			break;
		}

		if(n == (ApList.size() - 1))
		{
			cout << "No aperture found for InterpolatedEllipticalAperture at z = " << z << endl;
			for(size_t m = 0; m < ApList.size(); m++)
			{
				cout << "Entry " << m << " - z = " << ApList[m].s << "\t" << ApList[m].ap1
						<< endl;
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

	if(((x * x) / (ellipse_half_horizontal * ellipse_half_horizontal)) + ((y * y) / (ellipse_half_vertical
			* ellipse_half_vertical)) > 1)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

InterpolatedOctagonalAperture::InterpolatedOctagonalAperture(DataTable dt):
				InterpolatedAperture(dt)
{

}

InterpolatedOctagonalAperture::~InterpolatedOctagonalAperture()
{

}

Aperture* InterpolatedOctagonalAperture::GetInstance(DataTable dt)
{
	return new InterpolatedOctagonalAperture(dt);
}

string InterpolatedOctagonalAperture::GetType()
{
	return "OCTAGON-interpolated";
}

bool InterpolatedOctagonalAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	apStruct apFront;
	apStruct apBack;

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
	for(size_t n = 1; n < ApList.size(); n++)
	{
		if(ApList[n].s >= z)
		{
			apBack = ApList[n - 1];
			apFront = ApList[n];
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

	double hw = (g1 * z) + c1;
	double hh = (g2 * z) + c2;
	double angle1 = (g3 * z) + c3;
	double angle2 = (g4 * z) + c4;

	double tana1 = tan(angle1);
	double tana2 = tan(pi / 2 - angle2);

	double cc1 = (hh * tana2 - hw);
	double cc2 = hw * tana1;
	double cc3 = hh - cc2;

	//This is just taken from trrun.f90 in MAD-X. - credit to: 2015-Feb-20  18:42:26  ghislain: added octagon shape

	/*
		   !*** case of octagon: test outer rectangle (ap1,ap2) then test cut corner.
		   lost =  x .gt. ap1 .or. y .gt. ap2 .or. &
		     (ap2*tan(pi/2 - ap4) - ap1)*(y - ap1*tan(ap3)) - (ap2 - ap1*tan(ap3))*(x - ap1) .lt. zero
	 */
	double fabsx = fabs(x);
	double fabsy = fabs(y);

	x = fabsx;
	y = fabsy;

	if(x >= hw || y >= hh)
	{
		return false;
	}
	if(cc1 * (y - cc2) - cc3 * (x - hw) <= 0)
	{
		return false;
	}
	return true;
}

Aperture* InterpolatorFactory::GetInstance(DataTable dt)
{
	map<string, getInterpolator>::iterator itr = interpolatorTypes.find(dt.Get_s("APERTYPE",0));
	if(itr != interpolatorTypes.end())
	{
		return (*itr->second)(dt);
	}
	return NULL;
}

InterpolatorFactoryInitializer::InterpolatorFactoryInitializer()
{
	InterpolatorFactory::interpolatorTypes["RECTELLIPSE"] = &InterpolatedRectEllipseAperture::GetInstance;
	InterpolatorFactory::interpolatorTypes["CIRCLE"] = &InterpolatedCircularAperture::GetInstance;
	InterpolatorFactory::interpolatorTypes["RECTANGLE"] = &InterpolatedAperture::GetInstance;
	InterpolatorFactory::interpolatorTypes["ELLIPSE"] = &InterpolatedEllipticalAperture::GetInstance;
	InterpolatorFactory::interpolatorTypes["OCTAGON"] = &InterpolatedOctagonalAperture::GetInstance;
}

map<string, getInterpolator> InterpolatorFactory::interpolatorTypes;
InterpolatorFactoryInitializer InterpolatorFactoryInitializer::init;
