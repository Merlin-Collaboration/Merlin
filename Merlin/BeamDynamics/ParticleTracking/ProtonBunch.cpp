//#include <TH1.h>
#include <cmath>
#include "ProtonBunch.h"
#include "Collimators/CoulombScatter.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "Exception/MerlinException.h"
#include "AcceleratorModel/Aperture.h"
#include "AcceleratorModel/Apertures/TiltedAperture.hpp"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;


//extern TH1D* histt1;
//extern TH1D* histt2;

double ProtonBunch::GetParticleMass() const
{
	return ProtonMass;
}
double ProtonBunch::GetParticleMassMeV() const
{
	return ProtonMassMeV;
}

double ProtonBunch::GetParticleLifetime() const
{
	return 0;
}

bool ProtonBunch::IsStable() const
{
	//We assume protons do not decay :]
	return true;
}

int ProtonBunch::Scatter(PSvector& p,double x,double E0,const Aperture* ap)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference momentum
	// ap is the aperture

	//We simulate many physical processes that occur in bulk materials.
	// Ionization, multiple coulomb scattering,

	//distance travelled
	double z=0;

	// maximum allowed energy loss
	static const double MAXDP = 1.0 - 1.0e-7;

	//Check if we have a tilted aperture, aka a collimator jaw (with possible rotation).
	const TiltedAperture* tap= dynamic_cast<const TiltedAperture*> (ap);
	if(!tap) throw MerlinException("ScatterProton : No Tilted Aperture");

	//All is good, so lets extract scattering information from the aperture.
	double sig_pN_tot = tap->Material->sig_pN_tot;
	double sig_pN_in = tap->Material->sig_pN_in;
	double sig_R = tap->Material->sig_R;
	double dEdx = tap->Material->dEdx*tap->Material->rho/100;
	double rho = tap->Material->rho;
	double A = tap->Material->A;
	if (tap->Material->X0 == 0) throw MerlinException("Zero radiation length");
	double X0 = (tap->Material->X0*centimeter)/tap->Material->rho;

	double sig_pp_SD = 4.90e-3; // (barns) given in Chiara's thesis for 7TeV
	double sig_pp_el = 7.98e-3; // (barns) given in Chiara's thesis for 7TeV

	//Subtract a drift, this will soon be re-applied.
	p.x()-=x*p.xp();
	p.y()-=x*p.yp();

	// total mean free path (units meter)
	double lambda_tot = A * 1.e-6 / ((sig_pN_tot + sig_R) * barn* rho * Avogadro);

	// while there is remaining collimator
	while(x>0)
	{
	double t;
	double old_dp = p.dp();
	double delta_s = -lambda_tot * log (RandomNG::uniform (0, 1));
	double step_size = min(x, delta_s);
	double zstep=step_size*sqrt(1-p.xp()*p.xp()-p.yp()*p.yp());
	p.x()+=step_size*p.xp();
	p.y()+=step_size*p.yp();
	double E1 = E0 * (1+p.dp ());	// particle energy
	double thick = step_size / X0;	// material length in radiation lengths
	double dp = dEdx * step_size;	//(units of GeV )

	if(dp>MAXDP)
	{
		dp=MAXDP;
	}

	double E2 = E1 - dp;
	p.dp () =  (E2 - E0) / E0;

	// small-angle Coulomb scattering
	double Eav = (E1+E2)/2.0;
	double theta0 = 0.0136 * sqrt (thick) * (1.0 + 0.038 * log (thick)) / Eav;
	pair < double, double >s = CoulombScatter (step_size, theta0);
	p.x () += s.first;
	p.xp () += s.second;
	s = CoulombScatter (step_size, theta0);
	p.y () += s.first;
	p.yp () += s.second;
	if (E2 < (E0 / 100.0))
	{
		return 4;
	}

        if(tap->PointInside(p.x(),p.y(),z+=zstep))
        {
		tally[0]++; 
		p.x()+=p.xp()*x;
		p.y()+=p.yp()*x;
		return 0;
	}

	if (delta_s < x)
	{
		double E1 = E0 * (1+p.dp ());
		double r= RandomNG::uniform(0,1) * (sig_pN_tot + sig_R);
		double sig_pn_el =  1.6*pow(A,1/3.0)*sig_pp_el;
		double sig_pn_SD =  1.6*pow(A,1/3.0)*sig_pp_SD;
		double sig_pN_el = sig_pN_tot - sig_pN_in - sig_pn_el-sig_pn_SD;
		if ( (r -= sig_pN_el) < 0  )
		{
			// Elastic scatter pN (proton - Nucleus)
			tally[1]++;
			double b = 14.1*pow(A,2/3.0); //slope
			double mass=A*AtomicMassUnit;
			t = - log(RandomNG::uniform(0,1))/b;
                        //histt1->Fill(t);
			dp = t/(2*mass); // units of GeV
		}
		else if ( (r -= sig_pn_el) < 0  )
		{
                        tally[2]++;
			double b = 8.5 +1.086*log(114.59) ; // slope given on GeV units
			double tmax=0.011; //units of GeV^2
			while((t = - log(RandomNG::uniform(0,1))/b) > tmax);
                        //histt2->Fill(t);
			dp = t/(2*AtomicMassUnit); 
		}
		else if ( (r -= sig_pn_SD) < 0  )
		{
                        tally[3]++; 
			double b =6.5  ; //for Mx2 >2 GeV2
			double tmax= 0.1;
			double u=RandomNG::uniform(0,1);
			double Mx2;
			Mx2=pow((2.0*30.0)/(30.0-(u*28.0)),2);

			while (Mx2 < 2.0 || Mx2 > 30.0)
			{
				u=RandomNG::uniform(0,1); 
				Mx2=pow((2.0*30.0)/(30.0-(u*28.0)),2); 
			}
			
			t=- log(RandomNG::uniform(0,1))/b;

			while(t < tmax)
			{
				t=-log(RandomNG::uniform(0,1))/b;
			}
			dp = (t+Mx2-pow(AtomicMassUnit,2))/(2*AtomicMassUnit);
		}
		else if ( (r -= sig_R) < 0)
		{
                        tally[4]++;
			double tcut=0.000998;
			t=tcut/(1-RandomNG::uniform(0,1)); // generates 1/t squared distribution,
			double mass=A*AtomicMassUnit; // Nucleus mass
			dp = t/(2*mass); //units of GeV
		}
		else
		{
			tally[5]++; 
			p.dp()=-1.;

			//lose particle
			return 1;
		}

		double E2 = E1 - dp;
		p.dp()= (E2- E0)/E0;
		double theta=sqrt(t/(1+p.dp()))/E0;
		double phi=RandomNG::uniform(-pi,pi);
		p.xp()+=theta*cos(phi);
		p.yp()+=theta*sin(phi);
	}

	x -= step_size;
	} // end of while loop

	return 0;
} //End of ScatterProton
