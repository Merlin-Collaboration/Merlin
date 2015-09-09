#ifndef _h_ScatteringProcess
#define _h_ScatteringProcess 1

#include <iostream>
#include <cmath>

#include "merlin_config.h"

#include "BeamModel/PSvector.h"

#include "Collimators/Material.h"
#include "Collimators/DiffractiveScatter.h"
#include "Collimators/ElasticScatter.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

/*

Definition of the virtual ScatteringProcess class and
also several child classes derived from it
 
Created RJB 23 October 2012
Modified HR 07.09.2015

*/

namespace Collimation {
	
	struct CrossSections{
	private:
		double E0;
		double sig_pN_tot_ref;
		double sig_pN_inel_ref;
		double sig_R_ref;
		double sig_R;
		double sig_pp_tot;
		double sig_pp_el;
		double sig_pn_el;
		double sig_pp_sd;
		double sig_pn_sd;
		double sig_pN_tot;
		double sig_pN_inel;
		double sig_pN_el;
		double lambda_tot;
		double elastic_diff;
		double com_sqd;
		double density;
		double atomic_mass;
		double atomic_no;
		int scat_type;
		
	public:
		ParticleTracking::ppElasticScatter* ElasticScatter;
		ParticleTracking::ppDiffractiveScatter* DiffractiveScatter;
		
		void Set_E0(double a){E0 = a;}
		const double Get_E0(){return E0;}
		
		void Set_sig_pN_tot_ref(double a){sig_pN_tot_ref = a;}
		const double Get_sig_pN_tot_ref(){return sig_pN_tot_ref;}
		
		void Set_sig_pN_inel_ref(double a){ sig_pN_inel_ref = a;}
		const double Get_sig_pN_inel_ref(){return sig_pN_inel_ref;}
		
		void Set_sig_R_ref(double a){sig_R_ref = a;}
		const double Get_sig_R_ref(){return sig_R_ref;}
		
		void Set_sig_R(double a){sig_R = a;}
		const double Get_sig_R(){return sig_R;}
		
		void Set_sig_pp_tot(double a){sig_pp_tot = a;}
		const double Get_sig_pp_tot(){return sig_pp_tot;}
		
		void Set_sig_pp_el(double a){sig_pp_el = a;}
		const double Get_sig_pp_el(){return sig_pp_el;}
		
		void Set_sig_pn_el(double a){sig_pn_el = a;}
		const double Get_sig_pn_el(){return sig_pn_el;}
		
		void Set_sig_pp_sd(double a){sig_pp_sd = a;}
		const double Get_sig_pp_sd(){return sig_pp_sd;}
		
		void Set_sig_pn_sd(double a){sig_pn_sd = a;}
		const double Get_sig_pn_sd(){return sig_pn_sd;}
		
		void Set_sig_pN_tot(double a){sig_pN_tot = a;}
		const double Get_sig_pN_tot(){return sig_pN_tot;}
		
		void Set_sig_pN_inel(double a){sig_pN_inel = a;}
		const double Get_sig_pN_inel(){return sig_pN_inel;}
		
		void Set_sig_pN_el(double a){sig_pN_el = a;}
		const double Get_sig_pN_el(){return sig_pN_el;}
		
		void Set_lambda_tot(double a){lambda_tot = a;}
		const double Get_lambda_tot(){return lambda_tot;}
		
		void Set_elastic_diff(double a){elastic_diff = a;}
		const double Get_elastic_diff(){return elastic_diff;}
		
		void Set_com_sqd(double a){com_sqd = a;}
		const double Get_com_sqd(){return com_sqd;}
		
		void Set_density(double a){density = a;}
		const double Get_density(){return density;}
		
		void Set_atomic_mass(double a){atomic_mass = a;}
		const double Get_atomic_mass(){return atomic_mass;}
		
		void Set_atomic_no(double a){atomic_no = a;}
		const double Get_atomic_no(){return atomic_no;}
		
		void Set_scat_type(int a){scat_type = a;}
		const int Get_scat_type(){return scat_type;}

		
		//default constructor
		CrossSections(){
			Set_sig_pN_tot_ref(0.); 	
			Set_sig_pN_inel_ref(0.); 
			Set_sig_R_ref(0.);		
			Set_sig_R(0.);			
			Set_sig_pp_tot(0.);		
			Set_sig_pp_el(0.);		
			Set_sig_pn_el(0.);		
			Set_sig_pp_sd(0.);		
			Set_sig_pn_sd(0.);		
			Set_sig_pN_tot(0.);		
			Set_sig_pN_inel(0.);	
			Set_sig_pN_el(0.);		
			Set_scat_type(0.);			
			Set_lambda_tot(0.);		
			Set_elastic_diff(0.);		
			Set_com_sqd(0.); 			
			Set_density(0.);
			Set_atomic_mass(0.);
			Set_atomic_no(0.);
			ElasticScatter 		= NULL;
			DiffractiveScatter 	= NULL;
			Set_E0(0.);
		}
				
		//overloaded constructor
		CrossSections(Material* mat, double E, int scattertype){
			Set_sig_pN_tot_ref(mat->GetSixtrackTotalNucleusCrossSection());
			Set_sig_pN_inel_ref(mat->GetSixtrackInelasticNucleusCrossSection());
			Set_sig_R_ref(mat->GetSixtrackRutherfordCrossSection());
			Set_sig_R(0.);			
			Set_sig_pp_tot(0.);		
			Set_sig_pp_el(0.);		
			Set_sig_pn_el(0.);		
			Set_sig_pp_sd(0.);		
			Set_sig_pn_sd(0.);		
			Set_sig_pN_tot(0.);		
			Set_sig_pN_inel(0.);	
			Set_sig_pN_el(0.);		
			Set_scat_type(0.);			
			Set_lambda_tot(0.);		
			Set_elastic_diff(0.);		
			Set_com_sqd(0.); 	
			Set_density(mat->GetDensity()/1E3);
			Set_atomic_mass(mat->GetAtomicMass());
			Set_atomic_no(mat->GetAtomicNumber());
			Set_scat_type(scattertype);
			ElasticScatter 		= NULL;
			DiffractiveScatter 	= NULL;
			Set_E0(E);
			
			ConfigureCrossSections(Get_E0());		
			Set_lambda_tot(GetTotalMeanFreePath());		
			
			std::cout << "ScatteringProcess::CrossSections: dEdx = " << mat->GetSixtrackdEdx() << endl;
			std::cout << "ScatteringProcess::CrossSections: rho = " << mat->GetDensity()/1000.0 << endl;
			std::cout << "ScatteringProcess::CrossSections: A = " << mat->GetAtomicMass() << endl;
			std::cout << "ScatteringProcess::CrossSections: Z = " << mat->GetAtomicNumber() << endl;
			std::cout << "ScatteringProcess::CrossSections: X0 = " << mat->GetRadiationLengthInM() << endl;
			std::cout << "ScatteringProcess::CrossSections: b_N_ref = " << mat->GetSixtrackNuclearSlope() << endl;
		}
		
		inline bool operator==(const CrossSections& rhs){
			if( (this->sig_pN_tot_ref != rhs.sig_pN_tot_ref)	|| (this->scat_type != rhs.scat_type) || (this->lambda_tot != rhs.lambda_tot ))
			{ return 0;}
			else {return 1;}
		}	
		
		void ConfigureCrossSections(double E0){		
			if (scat_type == 4){ //Merlin
				Set_com_sqd( (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * E0) + (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV) );
				std::cout << "\n\nScatteringProcess::Configure: com_sqd = " << Get_com_sqd() << endl;
			}
			else { //SixTrack
				Set_com_sqd(2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * E0);	//ecmsq in SixTrack
				std::cout << "\n\nScatteringProcess::Configure: com_sqd = " << Get_com_sqd() << endl;
			}
			 
			double P_ref = 450.0 * PhysicalUnits::GeV;
			double pp_elastic_reference = 0.007;
			double pp_elastic_const = 0.04792;
			double pp_total_const = 0.05788;
			double free_nucleon_const = 1.618;
			double single_diffractive_const = 0.00068;
			double pp_tot_ref = 0.04;
			
			double free_nucleon_count =  free_nucleon_const * pow(atomic_mass,(1./3.));
			
			//SixTrack with Advanced Elastic Scattering					
			double sigma_pp_elasticEM = 0;		
			if((Get_scat_type() == 2) || (Get_scat_type() == 4))
			{
				ElasticScatter = new ParticleTracking::ppElasticScatter();
				ElasticScatter->SetTMin(1e-4);
				ElasticScatter->SetTMax(1.0);
				ElasticScatter->SetStepSize(1e-4);
				ElasticScatter->GenerateTDistribution(Get_E0());
			}
			
			//SixTrack with Advanced Single Diffractive Scattering			
			if((Get_scat_type() == 3) || (Get_scat_type() == 4))
			{
				//~ double sLHC = pow(114.6192348986088,2);
				double Mproton = 0.938272013;                            
				double Mpion = 0.1349766;
				double sLHC = (2*pow(Mproton,2)+2*Mproton*E0);
				double Mmin2 = pow(Mproton+Mpion,2);      
				double xi_th = Mmin2/sLHC; // (M_p + M_pion)^2/s 
				DiffractiveScatter = new ParticleTracking::ppDiffractiveScatter();			
				DiffractiveScatter->SetTMin(0.0001);
				DiffractiveScatter->SetTMax(4);
				DiffractiveScatter->SetTStepSize(1e-4);
				DiffractiveScatter->SetXiMin(xi_th);//Threshould at (M_proton + M_pion)^2/s
				DiffractiveScatter->SetXiMax(0.12);
				DiffractiveScatter->SetXiStepSize(1e-6);
				DiffractiveScatter->GenerateDistribution(Get_E0());
			}	

			//Merlin Scattering
			if(Get_scat_type() == 4){
				//total cross section calculation(ref. PDG) 
				const double Z_pp = 35.4548e-3;
				const double B_pp = 0.30810e-3;
				const double Y1_pp = 42.5323e-3;
				const double Y2_pp = 33.3433e-3;			
				const double eta1 = 0.45817;
				const double eta2 = 0.545;
				const double s0 = 28.998225*PhysicalUnits::GeV*PhysicalUnits::GeV;
				const double s1 =1*PhysicalUnits::GeV*PhysicalUnits::GeV;
				const double s = com_sqd;
				
				Set_sig_pp_tot(Z_pp + B_pp*pow(log (s/s0),2.0) + Y1_pp * pow(s1/s,eta1) -Y2_pp * pow(s1/s,eta2));
				
				Set_sig_pp_el(ElasticScatter->GetElasticCrossSectionN());
				sigma_pp_elasticEM = ElasticScatter->GetElasticCrossSection();
				Set_elastic_diff(sigma_pp_elasticEM - sig_pp_el);				
				
				Set_sig_pn_el(free_nucleon_count * sig_pp_el);
				
				Set_sig_pp_sd(DiffractiveScatter->GetDiffractiveCrossSection());
				
				Set_sig_pn_sd(free_nucleon_count * sig_pp_sd);
				
				Set_sig_pN_tot(sig_pN_tot_ref + free_nucleon_count * (sig_pp_tot - pp_tot_ref));
								
				Set_sig_pN_inel(sig_pN_inel_ref);
				
				Set_sig_pN_el(Get_sig_pN_tot() - Get_sig_pN_inel() - Get_sig_pn_el() - Get_sig_pn_sd());
				
				Set_sig_R(sig_R_ref);		
				
				std::cout << "\nScatteringProcess::Configure: sig_pN_tot_ref = " << Get_sig_pN_tot_ref() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_inel_ref = " << Get_sig_pN_inel_ref() << endl;
				std::cout << "ScatteringProcess::Configure: sig_R_ref = " << Get_sig_R_ref() << endl;				
				std::cout << "ScatteringProcess::Configure: sig_pn_el = " << Get_sig_pn_el() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pn_sd = " << Get_sig_pn_sd() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_inel = " << Get_sig_pN_inel() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_el = " << Get_sig_pN_el() << endl;
				std::cout << "ScatteringProcess::Configure: sig_R = " << Get_sig_R() << endl;
			}
			//SixTrack like scattering
			else{
				//pp total
				Set_sig_pp_tot(pp_tot_ref * pow((Get_E0() / P_ref), pp_total_const));
				
				//pp & pn elastic
				if(Get_scat_type() == 2){
					Set_sig_pp_el(ElasticScatter->GetElasticCrossSectionN());
					sigma_pp_elasticEM = ElasticScatter->GetElasticCrossSection();
					Set_elastic_diff(sigma_pp_elasticEM - sig_pp_el);					
				}
				else{
					Set_sig_pp_el(pp_elastic_reference * pow((Get_E0()/P_ref), pp_elastic_const));					
				}
				Set_sig_pn_el(free_nucleon_count * sig_pp_el);
				
				//pp & pn sd
				if(Get_scat_type() == 3){
					Set_sig_pp_sd(DiffractiveScatter->GetDiffractiveCrossSection());
				}
				else{
					Set_sig_pp_sd(single_diffractive_const * log(0.15 * Get_com_sqd() ));
				}
				Set_sig_pn_sd(free_nucleon_count * sig_pp_sd);
									
				//pN total
				Set_sig_pN_tot(sig_pN_tot_ref +  free_nucleon_count * (sig_pp_tot - pp_tot_ref));
				
				//pN inelastic
				Set_sig_pN_inel(sig_pN_inel_ref * sig_pN_tot / sig_pN_tot_ref);
				
				//pN elastic
				Set_sig_pN_el(sig_pN_tot - sig_pN_inel - sig_pn_el - sig_pn_sd);
				
				//Rutherford
				Set_sig_R(sig_R_ref);
				
				std::cout << "\nScatteringProcess::Configure: sig_pN_tot_ref = " << Get_sig_pN_tot_ref() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_inel_ref = " << Get_sig_pN_inel_ref() << endl;
				std::cout << "ScatteringProcess::Configure: sig_R_ref = " << Get_sig_R_ref() << endl;				
				std::cout << "ScatteringProcess::Configure: sig_pn_el = " << Get_sig_pn_el() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pn_sd = " << Get_sig_pn_sd() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_inel = " << Get_sig_pN_inel() << endl;
				std::cout << "ScatteringProcess::Configure: sig_pN_el = " << Get_sig_pN_el() << endl;
				std::cout << "ScatteringProcess::Configure: sig_R = " << Get_sig_R() << endl;
			}				
		
		}				
		
		double GetTotalMeanFreePath(){			
			//Merlin scattering
			if(Get_scat_type() == 4){
				Set_lambda_tot( (Get_atomic_mass() * 1E-6 / ( (Get_sig_pN_tot() + Get_atomic_no() * Get_elastic_diff()) * PhysicalUnits::barn * Get_density() * PhysicalConstants::Avogadro)) );
				//~ std::cout << "\n\tCrossSections::GetTotalMeanFreePath: Merlin config, lambda = " << Get_lambda_tot() << endl;
				return Get_lambda_tot();
			}
			//SixTrack + Advanced Elastic
			else if (Get_scat_type() == 2){
				Set_lambda_tot( (Get_atomic_mass() * 1E-6 / ( (Get_sig_pN_tot() + Get_sig_R()+ Get_elastic_diff()) * PhysicalUnits::barn * Get_density() * PhysicalConstants::Avogadro)) );
				// std::cout << "\n\tCrossSections::GetTotalMeanFreePath: ST + Adv. El config, lambda = " << Get_lambda_tot() << endl;
				return Get_lambda_tot();
			}
			//Sixtrack
			else{
				Set_lambda_tot( (Get_atomic_mass() * 1E-6 / ( (Get_sig_pN_tot() + Get_sig_R()) * PhysicalUnits::barn * Get_density() * PhysicalConstants::Avogadro)) );				
				//~ std::cout << "\tCrossSections::GetTotalMeanFreePath: SixTrack config, lambda = " << Get_lambda_tot() << endl;				
				return Get_lambda_tot();
			}	
	}
};

class ScatteringProcess {
public:
	double sigma; 					// Integrated cross section for this process

protected:
	double E0;						// Reference energy
	Material* mat; 					// Material of the collimator being hit	
	Collimation::CrossSections* cs;	// CrossSections object holding all configured cross sections
	double t;						// Momentum transfer
	
public:
	// The first function must be provided for all child classes, and probably the second as well
	virtual bool Scatter(PSvector& p, double E)=0;
	virtual void Configure(Material* matin, CrossSections* CSin){mat=matin; cs=CSin;}
	virtual std::string GetProcessType() const {return "ScatteringProcess";}
};

// Rutherford
class Rutherford:public ScatteringProcess
{
	double tmin;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "Rutherford";}
};

class SixTrackRutherford:public ScatteringProcess
{
	double tmin;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "SixTrackRutherford";}
};

// Elastic pn
class Elasticpn:public ScatteringProcess
{
	double b_pp; //slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "Elastic_pn";}
};

class SixTrackElasticpn:public ScatteringProcess
{
	double b_pp; //slope
public:	
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "SixTrackElasic_pn";}	
};

// Elastic pN
class ElasticpN:public ScatteringProcess
{
	double b_N; //slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "Elastic_pN";}
};

class SixTrackElasticpN:public ScatteringProcess
{
	double b_N; //slope
public:	
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "SixTrackElasic_pN";}	
};

// Single Diffractive
class SingleDiffractive:public ScatteringProcess
{
	double m_rec; //recoil mass
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "SingleDiffractive";}	
	
};

class SixTrackSingleDiffractive:public ScatteringProcess
{
	double m_rec; //recoil mass
	double dp;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "SixTrackSingleDiffractive";}	
	
};

// Inelastic
class Inelastic:public ScatteringProcess
{
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "Inelastic";}
};

} //end namespace Collimation

#endif
