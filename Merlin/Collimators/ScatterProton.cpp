#include <cmath>
#include <fstream>
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/utils.h"

#include "BeamModel/PSTypes.h"
#include "Collimators/ScatterProton.h"
#include "Random/RandomNG.h"

#include "AcceleratorModel/Apertures/TiltedAperture.hpp"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

// calculate random small-angle Coulomb scattering
pair < double, double > CoulombScatterp (double x, double theta0)
{
	// x - material lenght in meters
	// theta0 - RMS scattering angle (plane)
	// See particle data book section 27.3
	static const double root12 = sqrt (12.0);

	double z1 = RandomNG::normal (0, 1);
	double z2 = RandomNG::normal (0, 1);

	double theta_plane = z2 * theta0;
	double y_plane = z1 * x * theta0 / root12 + x * theta_plane / 2;
	return make_pair (y_plane, theta_plane);
}

int ScatterProton (PSvector& p, double x, double E0, const Aperture * tap)
{
	std::cout << "Scatter proton\t" << E0 << std::endl;

	//Check for possible badness
	if(tap == NULL)
	{
		std::cerr << "tap is NULL in ScatterProton.cpp - exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	if(tap->Material == NULL)
	{
		std::cerr << "Material is NULL in ScatterProton.cpp" << std::endl;
		exit(EXIT_FAILURE);
	}

	double sig_pN_tot 	= tap->Material->sig_pN_tot;
	double sig_pN_in 	= tap->Material->sig_pN_in;
	double sig_R 		= tap->Material->sig_R;
	double dEdx 		= tap->Material->dEdx;
	double rho 		= tap->Material->rho;
	double A 		= tap->Material->A;
	double X0 		= tap->Material->X0;
	double density 		= tap->Material->density;

	double sig_pp_SD = 4.90e-3; // (barns) given in Chiara's thesis for 7TeV
	double sig_pp_el = 7.98e-3; // (barns) given in Chiara's thesis for 7TeV
 
	if (isnan(X0) || X0 == 0)
	{
		cout << "X0 is zero in ScatterProton, this is very bad" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "X0:\t" << X0 << endl;

	// total mean free path - 1e3 factor to convert to kg.
	double lambda_tot = A * 1e-3 / ((sig_pN_tot + sig_R) * barn * density * Avogadro);

	// while there is remaining collimator
	while(x>0)
	{

	double delta_s = -lambda_tot * log (RandomNG::uniform (0, 1));
	cout << "length = " << x << "\tlambda_tot = " << lambda_tot << "\tdelta_s = " << delta_s << endl;
	cout << "lambda_tot:\t"<< lambda_tot << "\tdelta_s:\t" << delta_s << "\tx:\t" << x << endl;

	if (delta_s > x)
	{
		//"coulomb" scatter
		static const double MAXDP = 1.0 - 1.0e-7;	// maximum allowed energy loss
		double E1 = E0 * (1 + p.dp ());			// particle energy
		double t = x / X0;				// material length in radiation lengths
		double dp = dEdx * rho * x * 100;		//(units: MeV/g/cm2 *g/cm3 *100cm)

		if(dp > MAXDP)
		{
			dp = MAXDP;
		}

		// Adjust particle dp/p accordingly
		// p.dp() -= dp*(1.0+p.dp());

		// New particle (absolute) energy
		double E2 = E1 - (dp * MeV);
		p.dp () = (E2 - E0) / E0;

		// Compute the random angle scatter assuming that the particle has
		// an average energy of (E+E0)/2 in the collimator
		double Eav = (E1 + E2) / 2.0;
		cout << "Eav:\t" << Eav << endl;
		double theta0 = 13.6*MeV  * sqrt (t) * (1.0 + 0.038 * log (t)) / Eav;	// small-angle Coulomb scattering

		// cout <<" E2 and theta0 "<<E2<<" "<<theta0<<endl;     
		pair < double, double >s = CoulombScatterp (x, theta0);
		p.x () += s.first;
		p.xp () += s.second;

		s = CoulombScatterp (x, theta0);
		p.y () += s.first;
		p.yp () += s.second;

		// if particle has lost 99% of its energy, make it walk the plank
		if (E2 < (E0 / 100.0))
		{
			cout << "I've lost all my energy" << endl;
			return 4;
		}

		return 0;			//keep particle
	}

	else
	{
		double r = RandomNG::uniform(0,1) * (sig_pN_tot + sig_R);    
		double sig_pN_el = sig_pN_tot - sig_pN_in - 1.6*pow(A,1/3)*(sig_pp_SD + sig_pp_el); // = total - inelastic -quasielastic
		double sig_pn_el =  1.6*pow(A,1/3)*sig_pp_el;
		double sig_pn_SD =  1.6*pow(A,1/3)*sig_pp_SD;
		cout << "r = "<< r << "sig_in = " << sig_pN_in << "sig_pN_el= "<<sig_pN_el<<"sig_pn_el= "<<sig_pn_el<< "sig_pn_SD= "<<sig_pn_SD<<"sig_R= "<<sig_R<<endl;
		
		//In each case, subtract sigma from r, and check if it has gone negative
		//Elastic proton-nucleus scattering.
		if ( (r -= sig_pN_el) < 0  )
		{
			// Elastic scatter pN (proton - Nucleus)
			double b = 14.1 * pow(A,2/3); //slope    
			double mass=A*AtomicMassUnit; // AtomicMassUnit (MeV) is given in # "NumericalUtils/PhysicalConstants.h" 
			cout << "mass= "<< A << "\t" << mass << " E0 = " << E0 << endl; 
			double t = - log(RandomNG::uniform(0,1))/b; 
			p.dp()-=-(t+(pow(mass,2)+pow(AtomicMassUnit,2))/(2*AtomicMassUnit*E0));
			double theta = sqrt(t/(1+p.dp()))/E0;  // NEEDS CHECKING
			double phi=RandomNG::uniform(-pi,pi);
			p.x()+=x;
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*cos(phi);
			cout << "Elastic scatter off nucleus t = " << t << " angle= " << theta << endl;
			cout << "1" << endl;      
			//return 0; //keep particle
		}

		//Inelastic proton-nucleus scattering, the incident proton is destroyed and hence lost.
		else if ( (r -= sig_pN_in) < 0  )
		{
			// Inelastic scatter pN (proton - Nucleus)
			cout << "2" << endl;
			return 1; //loose particle
		}

		//elastic proton-nucleon scattering.
		else if ( (r -= sig_pn_el) < 0  )
		{
			// Elastic scatter pn (proton - nucleon)
			double b = 8.5 +1.086*log(114.59) ; // slope given on GeV units
			cout << " b= " << b <<endl;   
			double t;
			double tmax=0.011;
			while((t = - log(RandomNG::uniform(0,1))/b) > tmax);
			double phi=RandomNG::uniform(-pi,pi);
			p.dp()-=t/(2*AtomicMassUnit*E0);
			double theta = sqrt(t/(1+p.dp()))/E0;
			p.x()+=x;
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*cos(phi);
			cout << "Elastic scatter off nucleon t = " << t << " angle= " << theta << endl;
			cout << "3" << endl;
			//return 0; //keep particle
		}

		//Single diffractive scattering.
		else if ( (r -= sig_pn_SD) < 0  )
		{
			// Quasielastic "single difractive" scatter pn (proton - nucleon)
			double b = 6.5  ; //for Mx2 > 2 GeV^2
			double tmax= 0.1;
			double Mx2;
			double u = RandomNG::uniform(0,1);
			Mx2 = ( 2 * 5 ) / ( 5 - (u * 3));

			while (Mx2 < 2 || Mx2 > 5)
			{
				u = RandomNG::uniform(0,1);
				Mx2 = (2*5) / ( 5 - (u * 3));
			}// 2 <= Mx2 <= 5 as a function of 1/u

			double t = - log(RandomNG::uniform(0,1))/b;
			while(t > tmax)
			{
				t = -log(RandomNG::uniform(0,1))/b;
			}
			
			double phi=RandomNG::uniform(-pi,pi);
			p.dp()-=-(t+Mx2+pow(AtomicMassUnit,2))/(2*AtomicMassUnit*E0);

			double theta = sqrt(t/(1+p.dp()))/E0;
			p.x()+=x;
                     
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*cos(phi);
			cout << "Single diffractive t = " << t << " angle = " << theta << endl;
			cout << "4" << endl;
			//return 0; //keep particle
		}

		 // Rutherford Scattering
		else if ( (r -= sig_R) < 0  )
		{
			double tcut = 0.000998;
			double t = tcut/(1-RandomNG::uniform(0,1)); // generates 1/t squared distribution, we think
			double phi = RandomNG::uniform(-pi,pi);
			double mass = A*AtomicMassUnit; // Nucleus mass
			p.dp()-=t/(2*mass*E0);
			double theta=sqrt(t/(1+p.dp()))/E0;
			p.x()+=x;
                      
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*cos(phi);
			cout << "Rutherford scatter" << endl;
			cout << "5" << endl;
			//return 0; //keep the particle is the cattering angle is small
		}
	}

	x -= delta_s;
	}
	//cout << "Bad scatter - did not return correctly" << endl;
	cout << "reached end of colimator" << endl;
	return 0;
  
} //End of ScatterProton

/*
		double r = RandomNG::uniform (0, 1) * (sig_pN_tot + sig_R);
		double sig_pN_el = sig_pN_tot - sig_pN_in;
		// subtract sigma from r, and check if it has gone negative
		if ((r -= sig_pN_el) < 0)
		{
			cout << "1" << endl;
			return 1;		//loose particle
		}
		else if ((r -= sig_pN_in) < 0)
		{
			cout << "2" << endl;
			return 2;		//loose particle
		}
		else if ((r -= sig_R) < 0)
		{
			cout << "3" << endl;
			return 3;		//loose particle
		}
	}

	//Return 0 in case of fall through? - should not get here.
	cout << "Bad scatter - did not return correctly" << endl;
	return 0;
}
*/
