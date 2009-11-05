#include <cmath>
#include "BeamModel/PSTypes.h"
#include "Random/RandomNG.h"
#include <fstream>
#include "NumericalUtils/utils.h"

namespace {

	// calculate random small-angle Coulomb scattering
	pair<double,double> CoulombScatter(double x, double theta0)
	{
		// x - material lenght in meters
		// theta0 - RMS scattering angle (plane)
		// See particle data book section 27.3
		static const double root12 = sqrt(12.0);

		double z1 = RandomNG::normal(0,1);
		double z2 = RandomNG::normal(0,1);

		double theta_plane = z2*theta0;
		double y_plane = z1*x*theta0/root12+x*theta_plane/2;

		return make_pair(y_plane,theta_plane);
	}
};

void ScatterParticle(PSvector& p, double X0, double x, double E0)
{
	// compute the random energy loss (approximate formulas)
	// inputs:
	// p  - particle phase space vector   
	// X0 - material radiation length in meters
	// x  - material physical length in meters
	// E0 - reference energy in GeV

	static const double MAXDP=1.0-1.0e-7; // maximum allowed energy loss

	double E1=E0*(1+p.dp()); // particle energy

	double t  = x/X0; // material length in radiation lengths
	double t1 = t/log(2.0);
	double t2 = 0.5*((t1-1.0)/(t1+1.0));
	double ga = pow(t1*Gamma(t1),-1.0/t1);
	double gn = pow(RandomNG::uniform(0,1),1.0/t1)/ga;

	// relative energy loss (relative to E1)
	double dp = gn-t2*gn*gn+ga*(ga*(ga-1.0)+t2)*gn*gn*gn;
	if(dp>MAXDP)
		dp=MAXDP;

	// Adjust particle dp/p accordingly
	p.dp() -= dp*(1.0+p.dp());

	// New particle (absolute) energy
	double E2 = E0*(1.0+p.dp());

	// Compute the random angle scatter assuming that the particle has
	// an average energy of (E+E0)/2 in the collimator
	double Eav = (E1+E2)/2.0;
	double theta0 = 0.0136*sqrt(t)*(1.0+0.038*log(t))/Eav; // small-angle Coulomb scattering

	pair<double,double> s = CoulombScatter(x,theta0);
	p.x() += s.first;
	p.xp()+= s.second;

	s = CoulombScatter(x,theta0);
	p.y() += s.first;
	p.yp()+= s.second;
}

