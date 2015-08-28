#include <iostream>
#include <fstream>
#include <iomanip>

#include "Collimators/ScatteringProcess.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

#include "Random/RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;
using namespace Collimation;

/*

  File contains a useful general purpose routine, and then
  two routines (Configure and Scatter) for a whole lot of standard processes
  Note that amu is in MeV

*/

void ScatterStuff(PSvector& p, double t, double E0){                
             double E1 = (p.dp() + 1) * E0;        
             double theta = sqrt(t)/E1;
             double phi = RandomNG::uniform(-pi,pi);                         
             p.xp() += theta * cos(phi);
             p.yp() += theta * sin(phi);
}

void ScatterStuff(PSvector& p, double t, double m, double E0){ // scatter PSvector by t off target m 
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
//Not used?
/*
void ScatterStuff(PSvector& p, double t, double m, double mrec, double E0){ 
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
*/
void ScatterStuff(double dp, PSvector& p, double t, double E0){ 
// for diffractive process where delta is calculated in the scatter function
			double E1 = (p.dp() + 1) * E0;
			double E2 = E1 - dp/E0;
			p.dp() = (E2/E0)-1;
			double theta = sqrt(t)/E2;
			double phi = RandomNG::uniform(-pi,pi);                         
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
}

// Rutherford
void Rutherford::Configure(Material* matin, CrossSections* CSin){ 
	ScatteringProcess::Configure(matin, CSin);
	tmin = 0.9982E-3; // DeMolaize thesis page 29 [GeV^2]
	t = tmin/(1-RandomNG::uniform(0,1));	
	sigma = cs->Get_sig_R();
}

bool Rutherford::Scatter(PSvector& p, double E){
	//~ std::cout << "ScatteringProcess::Rutherford::Scatter p.x() = " << p.x() << std::endl;
	double TargetMass = amu*mat->GetAtomicMass();
	ScatterStuff(p, t, TargetMass, E);	
	return true;
}
// ST Rutherford
void SixTrackRutherford::Configure(Material* matin, CrossSections* CSin){ 
	ScatteringProcess::Configure(matin, CSin);
	tmin = 0.9982E-3; // DeMolaize thesis page 29 [GeV^2]
	t = tmin/(1-RandomNG::uniform(0,1));	
	sigma = cs->Get_sig_R();
}

bool SixTrackRutherford::Scatter(PSvector& p, double E){	
	//~ std::cout << "ScatteringProcess::SixTrackRutherford::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::SixTrackRutherford::Scatter" << endl;
	ScatterStuff(p, t, E);	
	return true;
}

// Elastic pn
void Elasticpn::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_el();
	t = cs->ElasticScatter->SelectT();
}	
bool Elasticpn::Scatter(PSvector& p, double E){
	//~ std::cout << "ScatteringProcess::Elasticpn::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::Elasticpn::Scatter" << endl;
	ScatterStuff(p, t, amu, E);
	return true;
}
// ST Elasticpn
void SixTrackElasticpn::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);	
	sigma = cs->Get_sig_pn_el();	
}	
bool SixTrackElasticpn::Scatter(PSvector& p, double E){	
	//~ std::cout << "ScatteringProcess::SixTrackElasticpn::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::SixTrackElasticpn::Scatter" << endl;
	double com_sqd = 2 * ProtonMassMeV * MeV * E;	//ecmsq in SixTrack
	b_pp = 8.5 + 1.086 * log(sqrt(com_sqd)) ; // slope given on GeV units
	t = -log(RandomNG::uniform(0,1))/b_pp;
	ScatterStuff(p, t, E);
	return true;
}

// Elastic pN
void ElasticpN::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pN_el();	
	double b_N_ref = matin->GetSixtrackNuclearSlope();	
	b_N = b_N_ref * (cs->Get_sig_pN_tot()/cs->Get_sig_pN_tot_ref());
	t = -log(RandomNG::uniform(0,1))/b_N;	
}	
bool ElasticpN::Scatter(PSvector& p, double E){	
	//~ std::cout << "ScatteringProcess::ElasticpN::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::ElasticpN::Scatter" << endl;
	double TargetMass = amu*mat->GetAtomicMass();
	ScatterStuff(p, t, TargetMass, E);		
	return true;
}

// ST Elastic pN
void SixTrackElasticpN::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pN_el();
	double b_N_ref = matin->GetSixtrackNuclearSlope();	
	b_N = b_N_ref * (cs->Get_sig_pN_tot()/cs->Get_sig_pN_tot_ref());	
	t = -log(RandomNG::uniform(0,1))/b_N;	
}	
bool SixTrackElasticpN::Scatter(PSvector& p, double E){
	//~ std::cout << "ScatteringProcess::SixTrackElasticpN::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::SixTrackElasticpN::Scatter" << endl;
	ScatterStuff(p, t, E);		
	return true;
}

// SD
void SingleDiffractive::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_sd();
	
	std::pair<double,double>TM = cs->DiffractiveScatter->Select();
	t = TM.first;
	m_rec = TM.second;
}	
bool SingleDiffractive::Scatter(PSvector& p, double E){
	//~ std::cout << "ScatteringProcess::SingleDiffractive::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::SingleDiffractive::Scatter" << endl;	
	double com_sqd = (2 * ProtonMassMeV * MeV * E) + 2 * pow( (ProtonMassMeV * MeV),2);
	double dp = m_rec * m_rec * E / com_sqd;
	ScatterStuff(dp, p, t, E);
	return true;
}

// ST SD
void SixTrackSingleDiffractive::Configure(Material* matin, CrossSections* CSin){
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_sd();	
	std::cout << "ScatteringProcess::SixTrackSingleDiffractive::Scatter sigma = " << sigma << std::endl;
}	
bool SixTrackSingleDiffractive::Scatter(PSvector& p, double E){	
	//~ std::cout << "ScatteringProcess::SixTrackSingleDiffractive::Scatter p.x() = " << p.x() << std::endl;
	//~ std::cout << "\nScatteringProcess::SixTrackSingleDiffractive::Scatter" << endl;	
	double com_sqd = 2 * ProtonMassMeV * MeV * E;	//ecmsq in SixTrack
	double b_pp = 8.5 + 1.086 * log(sqrt(com_sqd)) ; // slope given on GeV units
	double xm2 = exp(RandomNG::uniform(0,1)*log(0.15*com_sqd));
	double b = 0.0;
	if(xm2 < 2.0)
	{
		b = 2 * b_pp;
	}
	else if(2.0 <= xm2 && xm2 <= 5.0)
	{
		b = (106.0 - 17.0 * xm2 ) * b_pp / 26.0;	
	}
	else if(xm2 > 5.0)
	{
		b = 7.0 * b_pp / 12.0;
	}
	t =-log(RandomNG::uniform(0,1))/b;
	dp = xm2*E/com_sqd;
	ScatterStuff(dp, p, t, E);
	return true;
}

//Inelastic
void Inelastic::Configure(Material* matin, CrossSections* CSin){ 
    ScatteringProcess::Configure(matin, CSin);
    sigma= cs->Get_sig_pN_inel();
}

bool Inelastic::Scatter(PSvector& p, double E){		
		//~ std::cout << "ScatteringProcess::SixTrackSingleDiffractive::Scatter p.x() = " << p.x() << std::endl;
		return false;
	} // Particle is lost
