#ifndef _h_CrossSections
#define _h_CrossSections 1

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

namespace Collimation {
	
class CrossSections{
public:

	//default constructor
	CrossSections();
			
	//overloaded constructor
	CrossSections(Material* mat, double E, int scattertype);
	
	inline bool operator==(const CrossSections& rhs){
		if( (this->sig_pN_tot_ref != rhs.sig_pN_tot_ref)	|| (this->scat_type != rhs.scat_type) || (this->lambda_tot != rhs.lambda_tot ))
		{ return 0;}
		else {return 1;}
	}	
	
	void ConfigureCrossSections(double E0);
	
	double GetTotalMeanFreePath();
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
};
}

#endif
