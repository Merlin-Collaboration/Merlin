#ifndef CoulombScatter_h
#define CoulombScatter_h 1
#include <cmath>
#include "RandomNG.h"
namespace ParticleTracking
{

/**
* calculate random small-angle Coulomb scattering
*/
pair<double,double> CoulombScatter(double x, double theta0)
{
	/**
	* x - material length in meters
	* theta0 - RMS scattering angle (plane)
	* See particle data book section 27.3
	*/
	static const double root12 = sqrt(12.0);

	double z1 = RandomNG::normal(0,1);
	double z2 = RandomNG::normal(0,1);

	double theta_plane = z2*theta0;
	double y_plane = z1*x*theta0/root12+x*theta_plane/2;

	return make_pair(y_plane,theta_plane);
}
}
#endif
