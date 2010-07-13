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
/*
	ofstream prescat("Output/prescat.dat", ios::app);
	prescat << p.x() << " " << p.xp() << " " << p.y() << " " << p.yp() << " " << p.dp()<<endl;
	prescat.close();
*/
	const double L=x; // L is length, x is length remaining
	//std::cout << "Scatter proton\t" << E0 <<std:: endl;
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
	//std::cout << tap->sig_pN_tot << std::endl;
	double sig_pN_tot = tap->Material->sig_pN_tot;
	double sig_pN_in = tap->Material->sig_pN_in;
	double sig_R = tap->Material->sig_R;
	double dEdx = tap->Material->dEdx;
	double rho = tap->Material->rho; //density [g/cm^3]
	double A = tap->Material->A;
	double X0 = tap->Material->X0;
	double sig_pp_SD = 4.90e-3; // (barns) given in Chiara's thesis for 7TeV
	double sig_pp_el = 7.98e-3; // (barns) given in Chiara's thesis for 7TeV
 
	if (X0 == 0)
	{
		cout << "X0 is zero, this is very bad" << endl;
		exit(EXIT_FAILURE);
	}


	//cout << "X0:\t" << X0 << endl;
	double lambda_tot = A * 1.e-6 / ((sig_pN_tot + sig_R) * barn* rho * Avogadro);	// total mean free path (units meter)
	
	// while there is remaining collimator
	while(x>0){
	double old_dp = p.dp();
	double delta_s = -lambda_tot * log (RandomNG::uniform (0, 1));
	//cout << "length=" << x << " lambda_tot="<<lambda_tot << " delta_s="<<delta_s<<endl;   
	//cout << "lambda_tot:\t"<< lambda_tot << endl;
	//cout << "density:\t"<< rho<< endl;
	
	// either tracking to end of magnet, or step, which ever is smaller
	double step_size = min(x, delta_s);
	
	//do drift tracking
	p.x()+=step_size*p.xp();
	p.y()+=step_size*p.yp();
	//cout << "Coulom scater over p.x= " << p.x() <<endl;
	//"\tdelta_s:\t" << delta_s << "\tx:\t" << x << endl;
	
	//{
	     //"coulomb" scatter

	     	
		static const double MAXDP = 1.0 - 1.0e-7;	// maximum allowed energy loss
		double E1 = E0 * (1+p.dp ());	// particle energy
		// cout <<"E1= "<<E1<<endl;
		double t = step_size / X0;	// material length in radiation lengths

	//	cout << "x = " << x << "delta_s =" << delta_s <<endl;
		double dp = dEdx * step_size;	//(units of GeV )
		// cout <<"dp= "<<dp<<endl;

	//	cout << "ELOSS: MAXDP=" << MAXDP << " dp=" << dp << endl;

		//if(dp > MAXDP)
		//{
		//	dp = MAXDP;
		//}

		// Adjust particle dp/p accordingly
		// p.dp() -= dp*(1.0+p.dp());

		// New particle (absolute) energy
		double E2 = E1 - dp;
		p.dp () =  (E2 - E0) / E0;

		// Compute the random angle scatter assuming that the particle has
		// an average energy of (E+E0)/2 in the collimator
		double Eav = (E1 + E2)/2.0;                                                                                                                                  //cout << "Eav:\t" << Eav << endl;
		double theta0 = 0.0136 * sqrt (t) * (1.0 + 0.038 * log (t)) / Eav;	// small-angle Coulomb scattering
		//double theta0 = 13.6e6 * sqrt (t) * (1.0 + 0.038 * log (t)) / Eav;	// small-angle Coulomb scattering
		// cout <<" E2 and theta0 "<<E2<<" "<<theta0<<endl;     
		pair < double, double >s = CoulombScatterp (step_size, theta0);
		p.x () += s.first;
		p.xp () += s.second;

		s = CoulombScatterp (step_size, theta0);
		p.y () += s.first;
		p.yp () += s.second;

		// if particle has lost 99% of its energy, make it walk the plank
		if (E2 < (E0 / 100.0))
		{
			cout << "I've lost all my energy" << endl;
			return 4;
		}

		//return 0;			//keep particle
	//}
	//  cout << "fell through without returning" << endl;
	// return 99;

	if (delta_s < x)
	{
		double E1 = E0 * (1+p.dp ());
	//	cout << "energy of particle E1 = :\t"<< E1  <<endl;	
		double r= RandomNG::uniform(0,1) * (sig_pN_tot + sig_R);    
		double sig_pN_el = sig_pN_tot - sig_pN_in - 1.6*pow(A,1/3.0)*(sig_pp_SD + sig_pp_el); // = total - inelastic -quasielastic
		double sig_pn_el =  1.6*pow(A,1/3.0)*sig_pp_el;
		double sig_pn_SD =  1.6*pow(A,1/3.0)*sig_pp_SD;
          //	cout << "r = :\t"<< r << "sig_pN_in = :\t" << sig_pN_in << "sig_pN_el= :\t"<<sig_pN_el <<"sig_R= :\t"<<sig_R<<endl;
	//	cout << "A = :\t"<< A << "sig_pn_SD = :\t" << sig_pn_SD << "sig_pn_el= :\t"<<sig_pn_el<<endl;
	//	cout<< "density= :\t" <<rho <<endl;
			// subtract sigma from r, and check if it has gone negative
		if ( (r -= sig_pN_el) < 0  )
		{
			// Elastic scatter pN (proton - Nucleus)
			double b = 14.1*pow(A,2/3.0); //slope    
			double mass=A*AtomicMassUnit; // AtomicMassUnit (MeV) is given in # "NumericalUtils/PhysicalConstants.h" 
	//		cout << "mass= "<< A<< " "<< mass << "E0= " << E0<<endl; 
			double t = - log(RandomNG::uniform(0,1))/b; 
			double dp = t/(2*mass); // units of GeV
			double E2 = E1 - dp;
			p.dp()= (E2- E0)/E0;
	//		cout << "t = :\t"<< t << "dp = :\t" <<dp << "p.dp = :\t"<< p.dp() <<endl;
			double theta = sqrt(t/(1+p.dp()))/E0;  // NEEDS CHECKING
			double phi=RandomNG::uniform(-pi,pi);
			
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*sin(phi);
			
	//		cout<<"Elastic scatter of nucleus t = "<<t<<" angle= "<<theta<<endl;
	//		cout << "1" << endl;      
			//return 0; //keep particle
		}
		else if ( (r -= sig_pN_in) < 0  )
		{
			// Inelastic scatter pN (proton - Nucleus)
	//		cout << "2" << endl;
	//		cout <<"Inelastic " <<endl;
			return 1; //loose particle
		}
		else if ( (r -= sig_pn_el) < 0  )
		{
			// Elastic scatter pn (proton - nucleon)
			double b = 8.5 +1.086*log(114.59) ; // slope given on GeV units
	//		cout << " b= " << b <<endl;   
			double t;
			double tmax=0.011; //units of GeV^2
			while((t = - log(RandomNG::uniform(0,1))/b) > tmax); 
			double phi=RandomNG::uniform(-pi,pi);
			double dp = t/(2*AtomicMassUnit); 
			double E2 = E1 - dp;
			p.dp()= (E2- E0)/E0;
			double theta = sqrt(t/(1+p.dp()))/E0;  
	//	        cout << "t = :\t"<< t << "dp = :\t" <<dp << "p.dp = :\t"<< p.dp() <<endl;
			p.xp()+=theta*cos(phi);
			p.yp()+=theta*sin(phi);
			
	//		cout<<"Elastic scatter of nucleon t = "<<t<<" angle= "<<theta<<endl;
	//		cout << "3" << endl;
			//return 0; //keep particle
		}
	
		else if ( (r -= sig_pn_SD) < 0  )
		{
			// Quasielastic "single difractive" scatter pn (proton - nucleon)
			double b =6.5  ; //for Mx2 >2 GeV2
			double tmax= 0.1;
			double u=RandomNG::uniform(0,1);
			double Mx2;
			Mx2=pow((2.0*30.0)/(30.0-(u*28.0)),2);
			while (Mx2 < 2.0 || Mx2 > 30.0)
			{   u=RandomNG::uniform(0,1); 
			    Mx2=pow((2.0*30.0)/(30.0-(u*28.0)),2); 
			}
                        double t;
                        t=- log(RandomNG::uniform(0,1))/b;
			while(t < tmax)
			{   
				t=- log(RandomNG::uniform(0,1))/b;
			}
			
			double phi=RandomNG::uniform(-pi,pi);
			
			double dp = (t+Mx2-pow(AtomicMassUnit,2))/(2*AtomicMassUnit);
			double E2 = E1 -dp;
			p.dp()= (E2-E0)/E0;
			double theta = sqrt(t/(1+p.dp()))/E0;
	//	        cout << "t = :\t"<< t << "dp = :\t" <<dp << "p.dp = :\t"<< p.dp() <<endl;
                        p.xp()+=theta*cos(phi);
		        p.yp()+=theta*sin(phi);
			
	//		cout<<"Single diffractive t = "<<t<<" angle= "<<theta<<endl;
	//		cout << "4" << endl;
			//return 0; //keep particle
		}

		else if ( (r -= sig_R) < 0)  
		{
			 // Rutherford Scattering
			double tcut=0.000998;
			double t=tcut/(1-RandomNG::uniform(0,1)); // generates 1/t squared distribution, 
			double phi=RandomNG::uniform(-pi,pi);
			double mass=A*AtomicMassUnit; // Nucleus mass
			
			double dp = t/(2*mass); //units of GeV
			double E2 = E1 - dp;
			p.dp()= (E2- E0)/E0;
			double theta=sqrt(t/(1+p.dp()))/E0;
			
          //              cout << "t = :\t"<< t << "dp = :\t" <<dp << "p.dp = :\t"<< p.dp() <<endl;
                        p.xp()+=theta*cos(phi);
			p.yp()+=theta*sin(phi);
			
		//	cout << "Rutherford scatter"<<t<<" angle= "<<theta<<endl;
		//	cout << "5" << endl;
			//return 0; //keep the particle is the cattering angle is small
		}
	}

	x -= delta_s;
	if (p.dp() > old_dp){
		cerr << "ENERGY GAIN old:" << old_dp<< "  new"<<p.dp() << endl;
		exit(EXIT_FAILURE);
	}
	}
	//cout << "Bad scatter - did not return correctly" << endl;
	//cout << "reached end of colimator" << endl;
	/*
	ofstream postscat("Output/postscat.dat", ios::app);
	postscat << p.x() << " " << p.xp() << " " << p.y() << " " << p.yp() << " " << p.dp()<<endl;
	postscat.close();
	*/
	//merlin will presumably now track the particle through a drift, but we have already done this tracking, so subtract the value that merlin will now add
	p.x()-=L*p.xp();
	p.y()-=L*p.yp();
	
	//p.x() = 0;    // FIXME these are a test
	//p.xp() = 1e-6;// FIXME
	return 0;
  
} //End of ScatterProton

