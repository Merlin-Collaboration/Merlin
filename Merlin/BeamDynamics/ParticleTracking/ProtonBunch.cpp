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

int ProtonBunch::Scatter(PSvector& p, double x, double E0, const Aperture* ap)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference momentum
	// ap is the aperture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		TODO
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

      
	//Keep track of distance along the collimator for aperture checking (aperture could vary with z)
	double z = 0; 	// distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05;	// maximum allowed energy loss - 95%

	//Do a cast to check if we have a "TiltedAperture"
	const TiltedAperture* tap= dynamic_cast<const TiltedAperture*> (ap);
	//  if((++counter %  10000) ==0) cout<<" Event "<<counter<<endl;
	if(!tap)
	{
		throw MerlinException("ScatterProton : No Tilted Aperture");
	}

	const double sigma_pN_total_reference = tap->Material->sigma_pN_total;		//Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tap->Material->sigma_pN_inelastic;	//Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tap->Material->sigma_Rutherford;		//Material reference Rutherford scattering cross section
	//const double dEdx = tap->Material->dEdx*tap->Material->rho/10;			//dE/dx - (GeV m^-1)
	const double dEdx = tap->Material->dEdx;					//dE/dx - (GeV m^-1)
	const double rho = tap->Material->rho;						//density (g /cm^3)
	const double A = tap->Material->A;						//Atomic mass
	//const double Z = tap->Material->atomic_number;					//Atomic number
	//const double X0 = (tap->Material->X0*centimeter)/tap->Material->rho;
	const double X0 = tap->Material->X0;
	const double b_N_ref = tap->Material->b_N;

	//We have now read the material properties, now to scale these if required to the current energy scale etc
	double center_of_mass_squared = 2 * ProtonMassMeV * MeV * E0;	//ecmsq in SixTrack

	//pp cross-sections and parameters for energy dependence scaling
	const double p_reference = 450.0 * GeV;			//Reference energy at which the scattering data is based on.	(pref)
	const double pp_total_reference = 0.04;			//Proton total cross section at reference energy		(pptref)
	const double pp_elastic_reference = 0.007;		//Proton elastic cross section at reference energy		(pperef)
	const double single_diffractive_constant = 0.00068;	//Single Diffractive constant					(sdcoe)
	const double pp_total_constant = 0.05788;		//pp total cross section constant				(pptco)
	const double pp_elastic_constant = 0.04792;		//pp elastic scattering cross section power constant		(ppeco)
	const double free_nucleon_constant = 1.618;		//free nucleon constant						(freeco)
	const double t_low_cut = 0.0009982;			//Rutherford scattering cut scale (GeV^2)			(tlcut)
	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME
	//bnref -  b_N slope - nucleus elastic FIXME

	//Atom radius
	//double atomic_radius = 1.2e-15 * pow(A,(1.0/3.0);	// In m, remember elsewhere if using area, to convert to barns
	
	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A,1.0/3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.

	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy
	double sigma_pp_total = pp_total_reference * pow((E0 / p_reference),pp_total_constant);		//(pptot)
	
	//Next we also scale the elastic proton-nucleon cross section to the current energy
	//Remember this is for a single nucleon, when applying to a nucleus
	double sigma_pp_elastic = pp_elastic_reference * pow((E0 / p_reference),pp_elastic_constant);	// demolaize equation 1.21 / catalan 3.7
	//And here the elastic cross section is adjusted to the number of nucleons in this material
	double sigma_pn_elastic	= free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	double sigma_pp_SingleDiffractive = single_diffractive_constant * log(0.15 * center_of_mass_squared);
	//And again scale to the number of nucleons
	double sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	double sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	double sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);
	
	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference * sigma_pN_total / sigma_pN_total_reference;

	//Caluclate the full nucleus elastic contribution
	double sigma_pN_elastic	= sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive; 

	//Work on slopes next
	double b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)) ; // slope given on GeV units

	double b_N = b_N_ref * (sigma_pN_total/sigma_pN_total_reference);

	//Finally calculate the mean free path
	double lambda_tot = A * 1.e-6 / ((sigma_pN_total + sigma_Rutherford) * barn * rho * Avogadro);	// total mean free path (units meter)

	//Track back a drift
	p.x() -= x * p.xp();
	p.y() -= x * p.yp();
/*
	cout << "MFP:\t"  << lambda_tot << "\t" << sigma_pN_total_reference << endl;
	cout << "Total:\t" << sigma_pN_total << endl;
	cout << "pN el:\t" << sigma_pN_elastic << endl;
	cout << "pN inel:\t" << sigma_pN_inelastic << endl;
	cout << "pn el:\t" << sigma_pn_elastic << endl;
	cout << "pn SD:\t" << sigma_pn_SingleDiffractive << endl;
	cout << "pN Ruth:\t" << sigma_Rutherford << endl;
	
	cout << b_N << "\t" << free_nucleon_count << endl;
	abort();
*/
	while( x > 0 )
	{
		bool interacted;	//Interacted on this step or not - true/false?
		double t=0.0;		//Momentum transfer
		double delta_s = -lambda_tot * log (RandomNG::uniform (0, 1));
		double step_size;

        	interacted = ( x > delta_s );
	        step_size = interacted ? delta_s : x;
	
		//Do MCS + dE/dx
        	double zstep = step_size * sqrt( 1 - p.xp()*p.xp() - p.yp()*p.yp() );
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());	// particle energy
		double thick = step_size / X0;	// material length in radiation lengths
		
		double dp = dEdx * step_size;	//(units of GeV )

		
		double E2 = E1 - dp;
		p.dp () =  ((E1 - dp) - E0) / E0;
		double Eav = (E1+E2) / 2.0;

		/*
		*
		*	MULTIPLE COULOMB SCATTERING
		*
		*/

		double theta0 = 13.6*MeV * sqrt (thick) * (1.0 + 0.038 * log (thick)) / Eav;	// small-angle Coulomb scattering

		pair < double, double > s = CoulombScatter (step_size, theta0);
		p.x ()  += s.first;
		p.xp () += s.second;

		s = CoulombScatter (step_size, theta0);
		p.y ()  += s.first;
		p.yp () += s.second;

		if (E2 < (E0 / 100.0))
		{
			return 4;
		}

		//Check we are still inside the collimator
		//If not, the particle leaves
		if(tap->PointInside(p.x(),p.y(),z+=zstep))
		{
			tally[0]++; 
			p.x() += p.xp()*x;
			p.y() += p.yp()*x;
			return 0;
		}

//		interacted = false;
	        if (interacted)
		{
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0,1) * (sigma_pN_total + sigma_Rutherford);    

			//Point process interaction
			//Choose which scattering process to do

			/*
			*
			*	Elastic scatter pN (proton - Nucleus)
			*
			*/
			if ( (r -= sigma_pN_elastic) < 0  )
			{

				tally[1]++;
				t = -log(RandomNG::uniform(0,1))/b_N;
//				t = -log(0.5)/b_N;
				//histt1->Fill(t);
				double mass=A*AtomicMassUnit;
				dp = t/(2*mass); // units of GeV
//dp = 0.0;
			}

			/*
			*
			*	Elastic scatter pn (proton - nucleon)
			*
			*/
			else if ( (r -= sigma_pn_elastic) < 0  )
			{

				tally[2]++; 
				t = -log(RandomNG::uniform(0,1))/b_pp;
				//histt2->Fill(t);
				dp = t/(2*AtomicMassUnit);
//				t=0.0;
//				dp=0.0;
			}

			/*
			*
			*	Single Diffractive
			*
			*/
			else if ( (r -= sigma_pn_SingleDiffractive) < 0 )
			{
				tally[3]++; 

				double Mx2 = exp( RandomNG::uniform(0,1) * log(0.15 * center_of_mass_squared) );
				double b_sd = 1.0;
				double p = E1  * (1.0 - Mx2/center_of_mass_squared);
				dp = E1 - p;

				if ( Mx2 < 2.0 )
				{
					b_sd = 2.0 * b_pp;
				}
				else if (( Mx2 >= 2.0 ) && ( Mx2 <= 5.0 ))
				{
					b_sd = (106.0 - (17.0*Mx2)) *  b_pp / 36.0;
				}
				else if ( Mx2 > 5.0 )
				{
					b_sd = 7.0 * b_pp / 12.0;
				}
				t = -log(RandomNG::uniform(0,1))/b_sd;

//				t=0.0;
//				dp=0.0;
				//histt4->Fill(t);
				//histt3->Fill(Mx2);
			}

			/*
			*
			*	Rutherford coulomb scattering
			*
			*/
			else if ( (r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut=t_low_cut;			
				t=tcut/(1-RandomNG::uniform(0,1)); // generates 1/t squared distribution,
				double mass=A*AtomicMassUnit; // Nucleus mass
				dp = t/(2*mass); //units of GeV
				//ruth = 2.607e-4 * exp(-t * 0.8561e3 * emr^2 ) * (atomic_number / t) ^2

//				t=0.0;
//				dp=0.0;
			}

			/*
			*
			*	Inelastic interaction - no more proton :(
			*
			*/
			else
			{
				tally[5]++; 
				//p.dp()=-1.;
				//lose particle
				p.ct() = z;
				return 1;
			}
			double E3 = E2 - dp;
			p.dp() = (E3 - E0)/E0;

			double theta = sqrt(t)/E3;
			double phi = RandomNG::uniform(-pi,pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);

			//double phi = RandomNG::uniform(-1,1);
			//p.xp() += theta * phi;
			//p.yp() += theta * phi;

			/* 
			// SIXTRACK METHOD
			double teta = sqrt(t)/E3;
			double va,vb,va2,vb2,r2 = 2;

			while( r2 > 1.0)
			{
				va = (2.0 * RandomNG::uniform(0,1))-1.0; // -1 to +1
				vb = RandomNG::uniform(0,1);
				va2 = va*va;
				vb2 = vb*vb;
				r2 = va2 + vb2;
			}
			p.xp() += teta * (2.0*va*vb) / r2;
			p.yp() += teta * (va2 - vb2) / r2;
			*/
		}

	x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP)	//dp cut should be at 95%
	{
		cout << "95% P loss" << endl;
		//p.dp()=-1;
		return 1;
	} 

	return 0;

} //End of ScatterProton

//data (mname(i),i=1,nrmat)/ 'Be','Al','Cu','W','Pb','C','C2' /
// in Cs and CsRef,1st index: Cross-sections for processes
// 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
// 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs
