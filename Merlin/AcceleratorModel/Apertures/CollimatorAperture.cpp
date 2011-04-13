#include "AcceleratorModel/Apertures/CollimatorAperture.hpp"
#include "MADInterface/MADInterface.h"
#include "Random/RandomNG.h"

CollimatorAperture::CollimatorAperture(double w,double h, double t, material* m, bool err):RectangularAperture(w,h), alpha(t),errors(err),aperture_error(0.0)
{
	set_material(m);
}

//Checks if particle is in or outside a defined aperture
bool CollimatorAperture::PointInside(double x,double y,double z) const
{
	if(errors)
	{
		//Should use/adjust GetFullWidth
		x += aperture_error * ((pow(z,2)/jaw_length)-z);
		y += aperture_error * ((pow(z,2)/jaw_length)-z);
	}

	double x1 = (x * cos(alpha)) - (y * sin(alpha));
	double y1 = (x * sin(alpha)) + (y * cos(alpha));
	return fabs(x1) < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
}

void CollimatorAperture::EnableErrors(bool enable)
{
	errors = enable;
}

void CollimatorAperture::SetErrors(double max, double min)
{
	//if(errors)
	aperture_error = RandomNG::uniform(max,min);
	std::cout << "ApertureError:\t" << aperture_error << std::endl;
}

void CollimatorAperture::SetJawLength(double length)
{
	jaw_length = length;
}

bool CollimatorAperture::PointInside_offset(double x,double y,double z,const double xoff,const double yoff) const
{
	double x1 = (x-xoff) * cos(alpha);
	double y1 = (y-yoff) * sin(alpha);
	return fabs(x1)<GetFullWidth()/2 && fabs(y1)<GetFullHeight()/2;
}
