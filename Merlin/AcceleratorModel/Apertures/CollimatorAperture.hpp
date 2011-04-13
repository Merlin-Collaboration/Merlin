#ifndef _CollimatorAperture_hpp_
#define _CollimatorAperture_hpp_

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "Collimators/Material_Database.hpp"
#include "Collimators/Material.hpp"
#include <iostream>

class CollimatorAperture: public RectangularAperture
{

double alpha;
bool errors;
//double e1,e2,e3,e4,e5,e6;
double aperture_error;
double jaw_length;
public:

CollimatorAperture(double w,double h, double t, material* m, bool err=false);

void EnableErrors(bool);
void SetErrors(double, double);
void SetJawLength(double);

bool PointInside(double x,double y,double z) const;

bool PointInside_offset(double x, double y, double z, const double xoff, const double yoff) const;
};

#endif
