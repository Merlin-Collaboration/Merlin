#include "AcceleratorModel/Apertures/TiltedAperture.hpp"
#include "MADInterface/MADInterface.h"


TiltedAperture::TiltedAperture(double w,double h, double t, material* m):RectangularAperture(w,h), alpha(t)
{
	set_material( m);
};


//Checks if particle is in or outside a defined aperture
bool TiltedAperture::PointInside(double x,double y,double z) const
{
	double x1 = x * cos(alpha);
	double y1 = y * sin(alpha);
	return fabs(x1) < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
};

bool TiltedAperture::PointInside_offset(double x,double y,double z,const double xoff,const double yoff) const
{
	double x1 = (x-xoff) * cos(alpha);
	double y1 = (y-yoff) * sin(alpha);
	return fabs(x1)<GetFullWidth()/2 && fabs(y1)<GetFullHeight()/2;
};
