/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "merlin_config.h"

#include "PSvector.h"

#include "Material.h"
#include "CompositeMaterial.h"
#include "DiffractiveScatter.h"
#include "ElasticScatter.h"
#include "CrossSections.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;
using namespace Collimation;

//default constructor
CrossSections::CrossSections()
{
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
	ElasticScatter      = nullptr;
	DiffractiveScatter  = nullptr;
	Set_E0(0.);
}

//overloaded constructor
CrossSections::CrossSections(Material* mat, double E, int scattertype)
{
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
	Set_density(mat->GetDensity() / 1E3);
	Set_atomic_mass(mat->GetAtomicMass());
	Set_atomic_no(mat->GetAtomicNumber());
	Set_scat_type(scattertype);
	ElasticScatter      = nullptr;
	DiffractiveScatter  = nullptr;
	Set_E0(E);

	ConfigureCrossSections(Get_E0());
	Set_lambda_tot(GetTotalMeanFreePath());

	std::cout << "\nScatteringProcess::CrossSections: dEdx = " << mat->GetSixtrackdEdx() << endl;
	std::cout << "ScatteringProcess::CrossSections: rho = " << mat->GetDensity() / 1000.0 << endl;
	std::cout << "ScatteringProcess::CrossSections: A = " << mat->GetAtomicMass() << endl;
	std::cout << "ScatteringProcess::CrossSections: Z = " << mat->GetAtomicNumber() << endl;
	std::cout << "ScatteringProcess::CrossSections: X0 = " << mat->GetRadiationLengthInM() << endl;
	std::cout << "ScatteringProcess::CrossSections: b_N_ref = " << mat->GetSixtrackNuclearSlope() << endl;

	std::cout << "\nScatteringProcess::Configure: sig_pN_tot_ref = " << Get_sig_pN_tot_ref() << endl;
	std::cout << "ScatteringProcess::Configure: sig_pN_inel_ref = " << Get_sig_pN_inel_ref() << endl;
	std::cout << "ScatteringProcess::Configure: sig_R_ref = " << Get_sig_R_ref() << endl;

	std::cout << "\nScatteringProcess::Configure: sig_pN_tot = " << Get_sig_pN_tot() << endl;
	std::cout << "ScatteringProcess::Configure: sig_pN_inel = " << Get_sig_pN_inel() << endl;
	std::cout << "ScatteringProcess::Configure: sig_pN_el = " << Get_sig_pN_el() << endl;

	std::cout << "\nScatteringProcess::Configure: sig_pn_el = " << Get_sig_pn_el() << endl;
	std::cout << "ScatteringProcess::Configure: sig_pn_sd = " << Get_sig_pn_sd() << endl;
	std::cout << "ScatteringProcess::Configure: sig_R = " << Get_sig_R() << endl;

	std::cout << "\nScatteringProcess::CrossSections: Lambda_tot = " << GetTotalMeanFreePath() << endl;
}

CrossSections::~CrossSections()
{
	if(ElasticScatter)
	{
		delete ElasticScatter;
	}
	if(DiffractiveScatter)
	{
		delete DiffractiveScatter;
	}
}

void CrossSections::ConfigureCrossSections(double E0)
{
	if(scat_type == 4)  //Merlin
	{
		Set_com_sqd((2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * E0) + (2
			* PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * PhysicalConstants::ProtonMassMeV
			* PhysicalUnits::MeV));
		std::cout << "\n\nScatteringProcess::Configure: com_sqd = " << Get_com_sqd() << endl;
	}
	else   //SixTrack
	{
		Set_com_sqd(2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * E0);    //ecmsq in SixTrack
		std::cout << "\n\nScatteringProcess::Configure: com_sqd = " << Get_com_sqd() << endl;
	}

	double P_ref = 450.0 * PhysicalUnits::GeV;
	double pp_elastic_reference = 0.007;
	double pp_elastic_const = 0.04792;
	double pp_total_const = 0.05788;
	double free_nucleon_const = 1.618;
	double single_diffractive_const = 0.00068;
	double pp_tot_ref = 0.04;

	double free_nucleon_count =  free_nucleon_const * pow(atomic_mass, (1. / 3.));

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
		double sLHC = (2 * pow(Mproton, 2) + 2 * Mproton * E0);
		double Mmin2 = pow(Mproton + Mpion, 2);
		double xi_th = Mmin2 / sLHC; // (M_p + M_pion)^2/s
		DiffractiveScatter = new ParticleTracking::ppDiffractiveScatter();
		DiffractiveScatter->SetTMin(0.0001);
		DiffractiveScatter->SetTMax(4);
		DiffractiveScatter->SetTStepSize(1e-4);
		DiffractiveScatter->SetXiMin(xi_th); //Threshould at (M_proton + M_pion)^2/s
		DiffractiveScatter->SetXiMax(0.12);
		DiffractiveScatter->SetXiStepSize(1e-6);
		DiffractiveScatter->GenerateDistribution(Get_E0());
	}

	//Merlin Scattering
	if(Get_scat_type() == 4)
	{
		//total cross section calculation(ref. PDG)
		const double Z_pp = 35.4548e-3;
		const double B_pp = 0.30810e-3;
		const double Y1_pp = 42.5323e-3;
		const double Y2_pp = 33.3433e-3;
		const double eta1 = 0.45817;
		const double eta2 = 0.545;
		const double s0 = 28.998225 * PhysicalUnits::GeV * PhysicalUnits::GeV;
		const double s1 = 1 * PhysicalUnits::GeV * PhysicalUnits::GeV;
		const double s = com_sqd;

		Set_sig_pp_tot(Z_pp + B_pp * pow(log(s / s0), 2.0) + Y1_pp * pow(s1 / s, eta1) - Y2_pp * pow(s1 / s, eta2));

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

	}
	//SixTrack like scattering
	else
	{
		//pp total
		Set_sig_pp_tot(pp_tot_ref * pow((Get_E0() / P_ref), pp_total_const));

		//pp & pn elastic
		if(Get_scat_type() == 2)
		{
			Set_sig_pp_el(ElasticScatter->GetElasticCrossSectionN());
			sigma_pp_elasticEM = ElasticScatter->GetElasticCrossSection();
			Set_elastic_diff(sigma_pp_elasticEM - sig_pp_el);
		}
		else
		{
			Set_sig_pp_el(pp_elastic_reference * pow((Get_E0() / P_ref), pp_elastic_const));
		}
		Set_sig_pn_el(free_nucleon_count * sig_pp_el);

		//pp & pn sd
		if(Get_scat_type() == 3)
		{
			Set_sig_pp_sd(DiffractiveScatter->GetDiffractiveCrossSection());
		}
		else
		{
			Set_sig_pp_sd(single_diffractive_const * log(0.15 * Get_com_sqd()));
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
	}

}

double CrossSections::GetTotalMeanFreePath()
{
	//Merlin scattering
	if(Get_scat_type() == 4)
	{
		Set_lambda_tot((Get_atomic_mass() * 1E-6 / ((Get_sig_pN_tot() + Get_atomic_no() * Get_elastic_diff())
			* PhysicalUnits::barn * Get_density() * PhysicalConstants::Avogadro)));
		return Get_lambda_tot();
	}
	//SixTrack + Advanced Elastic
	else if(Get_scat_type() == 2)
	{
		Set_lambda_tot((Get_atomic_mass() * 1E-6 / ((Get_sig_pN_tot() + Get_sig_R() + Get_elastic_diff())
			* PhysicalUnits::barn * Get_density() * PhysicalConstants::Avogadro)));
		return Get_lambda_tot();
	}
	//Sixtrack
	else
	{
		Set_lambda_tot((Get_atomic_mass() * 1E-6 / ((Get_sig_pN_tot() + Get_sig_R()) * PhysicalUnits::barn
			* Get_density() * PhysicalConstants::Avogadro)));
		return Get_lambda_tot();
	}
}
