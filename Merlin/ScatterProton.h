#include "PSvector.h"
#include "MADInterface.h"
#include "utils.h"
#include "TiltedAperture.hpp"

pair<double,double> CoulombScatterp(double x, double theta0);
int ScatterProton(PSvector& p, double x, double E0,const  Aperture* tap);
//int ScatterProtonQ(PSvectorQ& p, double x, double E0,const  Aperture* tap);
