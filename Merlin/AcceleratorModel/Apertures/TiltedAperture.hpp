#ifndef _TiltedAperture_hpp_
#define _TiltedAperture_hpp_

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "Collimators/Material_Data.hpp"

class TiltedAperture: public RectangularAperture
{

double alpha;
int material;

public:
double sig_pN_tot;
double sig_pN_in;
double sig_R;
double dEdx;
double rho;
double sigma;
double A;
double X0;
//double tot_mean_free_path;
	
TiltedAperture(double w,double h, double t, int m=0):RectangularAperture(w,h), alpha(t), material(m)
{
	sig_pN_tot = coll_mat_data_sig_pN_tot[material];
	sig_pN_in = coll_mat_data_sig_pN_in[material];
	sig_R = coll_mat_data_sig_R[material];
	dEdx = coll_mat_data_dEdx[material];
	rho = coll_mat_data_rho[material];
	A = coll_mat_data_A[material];
	X0=coll_mat_data_X0[material];
	sigma = coll_mat_data_sigma[material];
	//tot_mean_free_path = 1 / sig_pN_tot;
};

virtual bool PointInside(double x,double y,double z) const;

virtual bool PointInside_offset(double x,double y,double z,const double xoff,const double yoff) const;
};

#endif
