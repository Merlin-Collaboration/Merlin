#include "BeamModel/PSvector.h"
#include "BeamModel/PSvectorQ.h"
#include "MADInterface/MADInterface.h"
#include "NumericalUtils/utils.h"
#include "AcceleratorModel/Apertures/TiltedAperture.hpp"

pair<double,double> CoulombScatterp(double x, double theta0);
int ScatterProton(PSvector& p, double x, double E0,const  Aperture* tap);
int ScatterProtonQ(PSvectorQ& p, double x, double E0,const  Aperture* tap);
