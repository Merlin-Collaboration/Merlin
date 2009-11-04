#include "BeamModel/PSvector.h"
#include "MADInterface/MADInterface.h"
#include "NumericalUtils/utils.h"
pair<double,double> CoulombScatterp(double x, double theta0);
int ScatterProton(PSvector& p, double X0, double x, double E0,const  TiltedAperture* tap);
