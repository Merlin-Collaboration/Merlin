#include <iostream>
#include <fstream>
#include <iomanip>

#include "Collimators/ScatteringModel.h"
#include "Collimators/ScatteringProcess.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

#include "Random/RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

/*

  File contains a useful general purpose routine, and then
  two routines (configure and scatter) for a whole lot of standard processes
  Note that amu is in MeV

*/

void scatterstuff(PSvector& p, double t, double m, double E0){ // scatter PSvector by t off target m 
// for delta in GeV (m[MeV] t[GeV^2])
             double delta = -t / (2*m);                          
             double E1 = (p.dp() + 1) * E0;
             double E2 = E1 - delta/E0;
             p.dp() = (E2/E0)-1;                        
             double theta = sqrt(t)/E2;
             double phi = RandomNG::uniform(-pi,pi);                         
             p.xp() += theta * cos(phi);
             p.yp() += theta * sin(phi);
}

void scatterstuff(PSvector& p, double t, double m, double mrec, double E0){ 
// for diffractive process where recoil mass mrec is larger than target mass m 
             double delta = ( -t + (pow(mrec,2) - pow(m,2)) ) / (m*2);
             double E1 = (p.dp() + 1) * E0;
			 double E2 = E1 - delta/E0;
             p.dp() = (E2/E0)-1;       		 
             double theta = sqrt(t)/E2;
             double phi = RandomNG::uniform(-pi,pi);                         
             p.xp() += theta * cos(phi);
             p.yp() += theta * sin(phi);
}

void scatterstuff(PSvector& p, double t, double m, double mrec, double E0, double dp){ 
// for diffractive process where delta is calculated in the scatter function
			double E1 = (p.dp() + 1) * E0;

			double E2 = E1 - dp/E0;
			p.dp() = (E2/E0)-1;
			//p.dp() = (E2 - E0) / E0;              
			//p.dp() -= dp;                         
			//cout << "\nscatterstuff3: E2 = " << E2 << " E1 = " << (p.dp() +1) * E0 << endl;
			double theta = sqrt(t)/E2;
			double phi = RandomNG::uniform(-pi,pi);                         
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
}

// Rutherford scattering functions
void Rutherford::Configure(Material* matin,double P0){ 
	ScatteringProcess::Configure(matin, P0);
	tmin=0.998E-3; // DeMolaize page 29 units [GeV]^2
	E0=P0;
	sigma = mat->GetSixtrackRutherfordCrossSection();
	}

bool Rutherford::scatter(PSvector& p, double E){
	double t=-tmin/(1-RandomNG::uniform(0,1));

	scatterstuff(p, t, mat->GetAtomicNumber()*amu, E);	
	return true;
	}


// Nucleon Elastic scattering functions, simple version
void NucleonElasticSimple::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    double sig_el_pp=0.001*7.0 * pow(P0/450.,.0479); //Lasheras Eq 3.7  converted to barns
    slope=8.5+1.086*log(sqrt(2*E0)); // Lasheras Eq 3.9
    double neff=1.6*pow(mat->GetAtomicNumber(),0.333); // Lasheras Eq 3.14
    sigma=neff*sig_el_pp;
}

bool NucleonElasticSimple::scatter(PSvector& p, double E){
	double t=-slope*log(RandomNG::uniform(0,1));	

	scatterstuff(p, t, ProtonMassGeV, E);	
	return true;
}


// Nucleus elastic scattering, sixtrack stylle

void NucleusElastic::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    sigma = mat->GetSixtrackElasticNucleusCrossSection();
    slope = 14.1*pow(mat->GetAtomicNumber(),0.65); // Lasheras Eq 3.15. Not quite right for mixtures
}

bool NucleusElastic::scatter(PSvector& p, double E){
	double t=-log(RandomNG::uniform(0,1))/slope;

	scatterstuff(p, t, ProtonMassGeV, mat->GetAtomicNumber()*amu, E);	
	return true; 
}


//Nucleon diffractive simple

void NucleonDiffractiveSimple::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    s = 2 * (ProtonMassGeV) * (E0 + (ProtonMassGeV));
    sigma=0.00068*log(s); // Lasheras 3.10
    Mx_lo2 = 3 * pow((ProtonMassGeV), 2);  // from Goulianos 
    Mx_hi2 =  Mx_lo2 + 0.15 * s;  // SD: Mx2_hi = Mp2 +0.15 *center_of_mass_squared
   }
  
bool NucleonDiffractiveSimple::scatter(PSvector& p, double E){
	double u = RandomNG::uniform(0,1);
	double Mx2 = Mx_lo2 * exp(u * log(Mx_hi2)); //Goulianos d2S/dMx2dt ~1/Mx2
	double b = 25.51 + 0.5 * log(s/Mx2); //Goulianos
	double t =-log(RandomNG::uniform(0,1))/b;

	scatterstuff(p, t, ProtonMassGeV, Mx2, E);
	
	return true;
} 

//Inelastic

void Inelastic::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    sigma= mat->GetSixtrackInelasticNucleusCrossSection(); 
}


//HR Sixtrack like scattering

void SixNucleonElastic::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    double sig_el_pp = 0.001 * 7.0 * pow((P0/450.),.0479); //Lasheras Eq 3.7  converted to barns
    double free_nucleon_count = 1.618 * pow((matin->GetAtomicNumber()),(1/3));
    double sigma_pn_el = free_nucleon_count * sig_el_pp;
	sigma = sigma_pn_el;  
}

bool SixNucleonElastic::scatter(PSvector& p, double E){
    double com = 2 * (ProtonMassGeV) * E0;
    //double com = 2 * ProtonMassMeV * E0;
    double b_pp = 8.5 + 1.086 * log( sqrt(com) );    
    double t = (-log(RandomNG::uniform(0,1)))/b_pp;
    
    scatterstuff(p, t, 1, mat->GetAtomicNumber()*amu, E);    
	return true; 
}

void SixSingleDiffractive::Configure(Material* matin, double P0){ 
    ScatteringProcess::Configure(matin,P0);
    double free_nucleon_count = 1.618 * pow((matin->GetAtomicNumber()),(1/3));
    double sd_const = 0.00068;
	double s = 2 * (ProtonMassGeV) * (E0 + ProtonMassGeV);
    comsqd = 2 * (ProtonMassGeV) * E0;
    double sig_pp_sd = sd_const * log( 0.15 * comsqd );
    double sig_pn_sd = free_nucleon_count * sig_pp_sd;
	sigma = sd_const * log(0.15 * comsqd);
    double b_pp = 8.5 + 1.086 * log( sqrt(comsqd) );
    
    xm2 = exp( RandomNG::uniform(0,1) * log(0.15 * comsqd) );
    b = 0;
    
    if(xm2 < 2.0){
        b = 2 * b_pp;
    }
    else if(xm2 <= 5.0 && 2.0 <= xm2){
        b = ( 106.0 - 17.0 * xm2 ) * b_pp / 26.0;
    }
    else if(xm2 > 5.0){
        b = 7 * b_pp / 12.0;
    } 
}

bool SixSingleDiffractive::scatter(PSvector& p, double E){
    double dp = xm2 * E0 / comsqd;
	double t = -log(RandomNG::uniform(0,1))/b;
	
	scatterstuff(p, t, ProtonMassGeV, xm2, E, dp);
	return true;
}



