#ifndef _TiltedAperture_hpp_
#define _TiltedAperture_hpp_

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "Collimators/Material_Database.hpp"
#include "Collimators/Material.hpp"
#include <iostream>

class TiltedAperture: public RectangularAperture
{

double alpha;

public:

TiltedAperture(double w,double h, double t, material* m);

virtual bool PointInside(double x,double y,double z) const;

virtual bool PointInside_offset(double x, double y, double z, const double xoff, const double yoff) const;
};

#endif
