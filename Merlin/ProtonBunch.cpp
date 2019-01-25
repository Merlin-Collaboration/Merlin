/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <fstream>
#include "ProtonBunch.h"
#include "CoulombScatter.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "MerlinException.h"
#include "Collimator.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

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

void ProtonBunch::EnableScatteringPhysics(scatMode chooseScat)
{
	switch(chooseScat)
	{
	case SixTrack:
		ScatteringPhysicsModel = 0;
		break;
	case SixTrackIoniz:
		ScatteringPhysicsModel = 1;
		break;
	case SixTrackElastic:
		ScatteringPhysicsModel = 2;
		break;
	case SixTrackSD:
		ScatteringPhysicsModel = 3;
		break;
	case Merlin:
		ScatteringPhysicsModel = 4;
		break;
	default:
		std::cout << "You need to select a scattering physics model" << std::endl;
	}
}

int ProtonBunch::Scatter(Particle& p, double x, const Collimator* col)
{
	int returnvalue;
	if(ScatteringPhysicsModel == 0)
	{
		returnvalue = ScatterSixtrack(p, x, col);
	}
	else if(ScatteringPhysicsModel == 1)
	{
		returnvalue = ScatterSixtrackAdvancedIonization(p, x, col);
	}
	else if(ScatteringPhysicsModel == 2)
	{
		returnvalue = ScatterSixtrackAdvancedElastic(p, x, col);
	}
	else if(ScatteringPhysicsModel == 3)
	{
		returnvalue = ScatterSixtrackAdvancedSingleDiffraction(p, x, col);
	}
	else if(ScatteringPhysicsModel == 4)
	{
		returnvalue = ScatterMerlin(p, x, col);
		MERLIN_PROFILE_END_TIMER("ProtonBunch::ScatterMerlin");
	}
	else
	{
		std::cerr << "Unknown scatter type: " << ScatteringPhysicsModel << std::endl;
		abort();
	}
	if(std::isnan(p.x()) || std::isnan(p.xp()) || std::isnan(p.y()) || std::isnan(p.yp()) || std::isnan(p.dp()) ||
		std::isnan(p.ct()))
	{
		std::cerr << "ProtonBunch::Scatter(): Particle has nan coordinate after scatter. ScatteringPhysicsModel="
				  << ScatteringPhysicsModel << std::endl;
		abort();
	}
	return returnvalue;
}
//End of ScatterProton

void ProtonBunch::ConfigureScatter(const Collimator* col)
{
	if(ScatteringPhysicsModel == 0)
	{
		ConfigureScatterSixtrack(col);
	}
	else if(ScatteringPhysicsModel == 1)
	{
		ConfigureScatterSixtrackAdvancedIonization(col);
	}
	else if(ScatteringPhysicsModel == 2)
	{
		ConfigureScatterSixtrackAdvancedElastic(col);
	}
	else if(ScatteringPhysicsModel == 3)
	{
		ConfigureScatterSixtrackAdvancedSingleDiffraction(col);
	}
	else if(ScatteringPhysicsModel == 4)
	{
		MERLIN_PROFILE_START_TIMER("ProtonBunch::ConfigureScatterMerlin");
		ConfigureScatterMerlin(col);
		MERLIN_PROFILE_END_TIMER("ProtonBunch::ConfigureScatterMerlin");
	}
}

void ProtonBunch::ConfigureScatterMerlin(const Collimator* col)
{

	//Do a cast to check if we have a "Collimator"
	const Collimator* tcol = dynamic_cast<const Collimator*>(col);
	if(!tcol)
	{
		std::cout << "ProtonBunch::ConfigureScatterMerlin() col is not Collimator" << std::endl;
		throw MerlinException("ScatterProton : No Collimator");
	}

	/**
	 * The proton momentum
	 */
	double P0 = GetReferenceMomentum();

	/**
	 * The proton total energy
	 */
	E0 = sqrt(P0 * P0 + pow(ProtonMassMeV * MeV, 2));

	/**
	 * The proton gamma factor
	 */
	double gamma = E0 / (ProtonMassMeV * MeV);

	/**
	 * The proton beta
	 */
	double beta = sqrt(1 - (1 / (gamma * gamma)));

	if(!GotElastic)
	{
		ElasticScatter = new ppElasticScatter();

		/**
		 * Generate the elastic differential cross section for the proton energy
		 */
		ElasticScatter->SetTMin(1e-4);
		ElasticScatter->SetTMax(1.0);
		ElasticScatter->SetStepSize(1e-4);
		ElasticScatter->GenerateTDistribution(E0);
		GotElastic = true;
	}

	if(!GotDiffractive)
	{
		double Mproton = 0.938272013;
		double Mpion = 0.1349766;
		double sLHC = (2 * pow(Mproton, 2) + 2 * Mproton * E0);

		double Mmin2 = pow(Mproton + Mpion, 2);
		double xi_th = Mmin2 / sLHC; // (M_p + M_pion)^2/s

		DiffractiveScatter = new ppDiffractiveScatter();

		/**
		 * Generate the Single Diffractive differential cross section for the proton energy
		 */
		DiffractiveScatter->SetTMin(0.0001);
		DiffractiveScatter->SetTMax(4);
		DiffractiveScatter->SetTStepSize(1e-4);
		DiffractiveScatter->SetXiMin(xi_th); //Threshold at (M_proton + M_pion)^2/s
		DiffractiveScatter->SetXiMax(0.12);
		DiffractiveScatter->SetXiStepSize(1e-6);

		DiffractiveScatter->GenerateDistribution(E0);
		GotDiffractive = true;
	}

	const double sigma_pN_total_reference = tcol->GetMaterial()->GetSixtrackTotalNucleusCrossSection();     //Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tcol->GetMaterial()->GetSixtrackInelasticNucleusCrossSection(); //Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tcol->GetMaterial()->GetSixtrackRutherfordCrossSection();     //Material reference Rutherford scattering cross section
	//const double dEdx = tcol->Material->dEdx*tcol->Material->rho/10;			//dE/dx - (GeV m^-1)
	//dEdx = tcol->GetMaterial()->dEdx;					//dE/dx - (GeV m^-1)
	rho = tcol->GetMaterial()->GetDensity() / 1000.0;                 //density (g /cm^3)
	A = tcol->GetMaterial()->GetAtomicMass();               //Atomic mass
	Z = tcol->GetMaterial()->GetAtomicNumber();             //Atomic number
	//const double X0 = (tcol->Material->X0*centimeter)/tcol->Material->rho;
	name = tcol->GetMaterial()->GetName();
	X0 = tcol->GetMaterial()->GetRadiationLengthInM();
	I = tcol->GetMaterial()->GetMeanExcitationEnergy() / eV;
	const double ElectronDensity = tcol->GetMaterial()->GetElectronDensity();       //N_e / m^3
	const double PlasmaEnergy = tcol->GetMaterial()->GetPlasmaEnergy() / eV;
	const double b_N_ref = tcol->GetMaterial()->GetSixtrackNuclearSlope();

	/**
	 * We have now read the material properties, now to scale these if required to the current energy scale etc
	 */
	center_of_mass_squared = (2 * ProtonMassMeV * MeV * E0) + (2 * ProtonMassMeV * MeV * ProtonMassMeV * MeV);

	//pp cross-sections and parameters for energy dependence scaling
	const double pp_total_reference = 0.04;         //Proton total cross section at reference energy		(pptref)
	const double free_nucleon_constant = 1.618;     //free nucleon constant						(freeco)
	t_low_cut = 0.0009982;                  //Rutherford scattering cut scale (GeV^2)			(tlcut)
	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME

	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A, 1.0 / 3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.

	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy

	//total cross section calculation(ref. PDG)
	const double Z_pp = 35.4548e-3;
	const double B_pp = 0.30810e-3;
	const double Y1_pp = 42.5323e-3;
	const double Y2_pp = 33.3433e-3;

	const double Z_pn = 35.8016e-3;
	const double B_pn = 0.30810e-3;
	const double Y1_pn = 40.15e-3;
	const double Y2_pn = 30.0e-3;

	const double eta1 = 0.45817;
	const double eta2 = 0.545;
	const double s0 = 28.998225 * GeV * GeV;
	const double s1 = 1 * GeV * GeV;
	const double s = center_of_mass_squared;

	//* Proton-proton total cross section
	double sigma_pp_total = Z_pp + B_pp * pow(log(s / s0), 2.0) + Y1_pp * pow(s1 / s, eta1) - Y2_pp * pow(s1 / s, eta2);

	//* Proton-neutron total cross section

	double sigma_pn_total = Z_pn + B_pn * pow(log(s / s0), 2.0) + Y1_pn * pow(s1 / s, eta1) - Y2_pn * pow(s1 / s, eta2);

	/**
	 * The total elastic cross section
	 */
	const double sigma_pp_elastic = ElasticScatter->GetElasticCrossSectionN();
	const double sigma_pp_elasticEM = ElasticScatter->GetElasticCrossSection();
	double ElasticDifference = sigma_pp_elasticEM - sigma_pp_elastic;

	//And here the elastic cross section is adjusted to the number of nucleons in this material
	sigma_pn_elastic = free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	double sigma_pp_SingleDiffractive = DiffractiveScatter->GetDiffractiveCrossSection();

	//And again scale to the number of nucleons
	sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);

	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference;

	//Caluclate the full nucleus elastic contribution
	sigma_pN_elastic = sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive;

	//Work on slopes next
	b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)); // slope given on GeV units

	b_N = b_N_ref * (sigma_pN_total / sigma_pN_total_reference);

	tmax = (2 * ElectronMassMeV * beta * beta * gamma * gamma) / (1 + (2 * gamma * (ElectronMassMeV / ProtonMassMeV))
		+ pow((ElectronMassMeV / ProtonMassMeV), 2)) * MeV;

	static const double xi1 = 2.0 * pi * pow(ElectronRadius, 2) * ElectronMass * pow(SpeedOfLight, 2);
	xi0 = xi1 * ElectronDensity;

	C = 1 + 2 * log(I / PlasmaEnergy);

	if((I / eV) < 100)
	{
		if(C <= 3.681)
		{
			C0 = 0.2;
			C1 = 2.0;
		}
		else
		{
			C0 = 0.326 * C - 1.0;
			C1 = 2.0;
		}
	}
	else    //I >= 100eV
	{
		if(C <= 5.215)
		{
			C0 = 0.2;
			C1 = 3.0;
		}
		else
		{
			C0 = 0.326 * C - 1.5;
			C1 = 3.0;
		}
	}

	//Finally calculate the mean free path
	lambda_tot = A * 1.e-6 / ((sigma_pN_total + Z * ElasticDifference) * barn * rho * Avogadro);  // total mean free path (units meter)

	SetScatterConfigured(true);
	bool output_scattering_details = false;
	if(output_scattering_details)
	{
		std::cout << "Merlin Scattering Configuration" << std::endl;
		std::cout << "..............................." << std::endl;
		std::cout << "Material : " << name << std::endl;
		std::cout << "Look at the material database reference:" << std::endl;
		std::cout << "Material reference proton-Nucleus total scattering cross section = "
				  << sigma_pN_total_reference << std::endl;
		std::cout << "Material reference proton-Nucleus inelastic scattering cross section = "
				  << sigma_pN_inelastic_reference << std::endl;
		std::cout << "Material reference Rutherford scattering cross section = " << sigma_Rutherford_reference
				  << std::endl;
		std::cout << "Density rho = " << rho << std::endl;                      //density (g /cm^3)
		std::cout << "Mass Number A = " << A << std::endl;
		std::cout << "Atomic number Z = " << Z << std::endl;
		std::cout << "X0 = " << X0 << std::endl;
		std::cout << "Nuclear reference slope b_N_ref = " << b_N_ref << std::endl;
		std::cout << "Now we need to scale the material properties to the required energy:" << std::endl;
		std::cout << "Center of mass energy squared = " << center_of_mass_squared << std::endl;
		std::cout << "Proton total cross section at reference energy (pptref) = " << pp_total_reference << std::endl;
		std::cout << "free nucleon constant = " << free_nucleon_constant << std::endl;
		std::cout << "Rutherford scattering cut scale (GeV^2) = " << t_low_cut << std::endl;
		std::cout << "First we calculate the number of nucleons available for scattering: " << std::endl;
		std::cout << "Number of free nucleons available for scattering= " << free_nucleon_count << std::endl;
		std::cout << "Scaling to the number of nucleons: " << std::endl;
		std::cout << "Scaled total pp cross section = " << sigma_pp_total << std::endl;
		std::cout << "Scaled total pn cross section = " << sigma_pn_total << std::endl;
		std::cout << "Total Elastic cross section = " << sigma_pp_elastic << std::endl;
		std::cout << "Scaled total Elastic cross section = " << sigma_pn_elastic << std::endl;
		std::cout << "Scaled Single diffractive cross section = " << sigma_pp_SingleDiffractive << std::endl;
		std::cout << "Scaled Single diffractive cross section Nucleus tot = " << sigma_pn_SingleDiffractive
				  << std::endl;
		std::cout << "sigma Rutherford(bug : not scaled)" << sigma_Rutherford << std::endl;
		std::cout << "Scaled total without Rutherford inclusion = " << sigma_pN_total << std::endl;
		std::cout << "Scaled Total inelastic = " << sigma_pN_inelastic << std::endl;
		std::cout << "Scaled full nucleus elastic contribution = " << sigma_pN_elastic << std::endl;
		std::cout << "b_pp(GeV units) = " << b_pp << std::endl;
		std::cout << "b_N (GeV units)= " << b_N << std::endl;
		std::cout << "Total mean free path (units meter) = " << lambda_tot << std::endl;
		std::cout << "..............................." << std::endl;
	}
}

int ProtonBunch::ScatterMerlin(PSvector& p, double x, const Collimator* col)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference energy
	// col is the colerture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		OK
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

	if(!ScatterConfigured)
	{
		ConfigureScatter(col);
	}
	MERLIN_PROFILE_START_TIMER("ProtonBunch::ScatterMerlin");
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);

	//Keep track of distance along the collimator for colerture checking (colerture could vary with z)
	double z = 0;   // distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05; // maximum allowed energy loss - 95%

	double TargetMass = A * AtomicMassUnit; // Nucleus mass

	while(x > 0)
	{
		bool interacted;    //Interacted on this step or not - true/false?
		double t = 0.0;       //Momentum transfer
		double delta_s = -lambda_tot * log(RandomNG::uniform(0, 1));
		double step_size;

		interacted = (x > delta_s);
		step_size = interacted ? delta_s : x;

		//Do MCS + dE/dx
		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());  // particle energy
		double thick = step_size / X0;  // material length in radiation lengths

		//Start dEdx
		double gamma = E1 / (ProtonMassMeV * MeV);
		double beta = sqrt(1 - (1 / (gamma * gamma)));

		double land = RandomNG::landau();

		double xi = (xi0 * step_size / (beta * beta)) / ElectronCharge * (eV / MeV);

		//Density correction
		double ddx = log10(beta * gamma);
		if(ddx > C1)
		{
			delta = 4.606 * ddx - C;
		}
		else if(ddx >= C0 && ddx <= C1)
		{
			double m = 3.0;
			double xa = C / 4.606;
			double a = 4.606 * (xa - C0) / pow((C1 - C0), m);
			delta = 4.606 * ddx - C + a * pow((C1 - ddx), m);
		}
		else
		{
			delta = 0.0;
		}

		double tcut = 2.0 * MeV;
		tcut = tmax;

		//Mott Correction
		double G = pi * FineStructureConstant * beta / 2.0;
		double q = (2 * (tmax / MeV) * (ElectronMassMeV)) / (pow((0.843 / MeV), 2));
		double S = log(1 + q);

		double L1 = 0.0;

		double yL2 = FineStructureConstant / beta;

		double L2sum = 1.202001688211;  //Sequence limit calculated with mathematica
		double L2 = -yL2 * yL2 * L2sum;

		double F = G - S + 2 * (L1 + L2);
		double deltaE = xi * (log(2 * ElectronMassMeV * beta * beta * gamma * gamma * (tcut / MeV) / pow(I / MeV, 2))
			- (beta * beta) * (1 + ((tcut / MeV) / (tmax / MeV))) - delta + F - 1.0 - euler);
		deltaE = xi * (log(2 * ElectronMassMeV * beta * beta * gamma * gamma * xi / pow(I / MeV, 2)) - (beta * beta)
			- delta + F + 0.20);

		if(xi > (2 * ElectronMassMeV * beta * beta * gamma * gamma))
		{
			cout << "de,xi " << deltaE << "\t" << xi << endl;
		}

		double dp = ((xi * land) - deltaE) * MeV;
		//End dEdx

		double E2 = E1 - dp;

		if(E2 <= 1.0)
		{
			p.ct() = z;
			return 1;
		}

		p.dp() = ((E1 - dp) - E0) / E0;
		double Eav = (E1 + E2) / 2.0;

		/*
		 *
		 *	MULTIPLE COULOMB SCATTERING
		 *
		 */
		double theta0 = 13.6 * MeV * sqrt(thick) * (1.0 + 0.038 * log(thick)) / Eav;    // small-angle Coulomb scattering

		/**
		 * Scatter in x
		 */
		pair<double, double> s = CoulombScatter(step_size, theta0);
		p.x() += s.first;
		p.xp() += s.second;

		/**
		 * Scatter in y
		 */
		s = CoulombScatter(step_size, theta0);
		p.y() += s.first;
		p.yp() += s.second;

		if(E2 < (E0 / 100.0))
		{
			return 4;
		}

		/**
		 * Check we are still inside the collimator
		 * If not, the particle leaves
		 */
		if(tcol->GetAperture()->CheckWithinApertureBoundaries(p.x(), p.y(), z += zstep))
		{
			tally[0]++;
			p.x() += p.xp() * x;
			p.y() += p.yp() * x;
			return 0;
		}

		/**
		 * Point process nuclear interactions
		 */
		if(interacted)
		{
			if(tcol->GetMaterial()->IsMixture())
			{
				tcol->GetMaterial()->SelectRandomMaterial();
			}
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0, 1) * (sigma_pN_total + sigma_Rutherford);

			//Choose which scattering process to do
			/**
			 *	Elastic scatter pN (proton - Nucleus)
			 */
			if((r -= sigma_pN_elastic) < 0)
			{
				tally[1]++;
				t = -log(RandomNG::uniform(0, 1)) / b_N;
				dp = t / (2 * TargetMass); // units of GeV
				p.type() = 1;
			}

			/**
			 *	Elastic scatter pn (proton - nucleon)
			 */
			else if((r -= sigma_pn_elastic) < 0)
			{
				tally[2]++;
				t = ElasticScatter->SelectT();
				dp = t / (2 * AtomicMassUnit);
				p.type() = 1;
			}

			/**
			 *	Single Diffractive
			 */
			else if((r -= sigma_pn_SingleDiffractive) < 0)
			{
				tally[3]++;

				std::pair<double, double> TM = DiffractiveScatter->Select();
				t = TM.first;
				double mrec = TM.second;

				dp = mrec * mrec * E1 / center_of_mass_squared;
				p.type() = 2;
				p.sd() = 1;
			}

			/**
			 *	Rutherford coulomb scattering
			 */
			else if((r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut = t_low_cut;
				t = tcut / (1 - RandomNG::uniform(0, 1)); // generates 1/t squared distribution,
				dp = t / (2 * TargetMass); //units of GeV
				p.type() = 3;
			}

			/**
			 *	Inelastic interaction - no more proton :(
			 */
			else
			{
				//lose particle
				tally[5]++;

				/**
				 * We use the p.ct() coordinate to set the exact position of
				 * the proton loss within the collimator - used for output.
				 */
				p.ct() = z;
				return 1;
			}

			double E3 = E2 - dp;

			if(E3 <= 0.10)
			{
				p.ct() = z;
				return 1;
			}

			p.dp() = (E3 - E0) / E0;

			double theta = sqrt(t) / E3;
			double phi = RandomNG::uniform(-pi, pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
		}

		x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP || p.dp() < -1.0)    //dp cut should be at 95%
	{
		p.ct() = z;
		return 1;
	}

	return 0;
}

/****************************************************************
 *
 *
 *	"Sixtrack/K2 style" scattering physics
 *
 *
 ****************************************************************/

void ProtonBunch::ConfigureScatterSixtrack(const Collimator* col)
{
	//Do a cast to check if we have a "Collimator"
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);
	if(!tcol)
	{
		throw MerlinException("ScatterProton : No Collimator colerture");
	}

	double P0 = GetReferenceMomentum();
	E0 = sqrt(P0 * P0 + pow(ProtonMassMeV * MeV, 2));

	const double sigma_pN_total_reference = tcol->GetMaterial()->GetSixtrackTotalNucleusCrossSection();      //Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tcol->GetMaterial()->GetSixtrackInelasticNucleusCrossSection();  //Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tcol->GetMaterial()->GetSixtrackRutherfordCrossSection();      //Material reference Rutherford scattering cross section
	dEdx = tcol->GetMaterial()->GetSixtrackdEdx();                   //dE/dx - (GeV m^-1)
	rho = tcol->GetMaterial()->GetDensity() / 1000.0;
	A = tcol->GetMaterial()->GetAtomicMass();                        //Atomic mass
	Z = tcol->GetMaterial()->GetAtomicNumber();              //Atomic number
	X0 = tcol->GetMaterial()->GetRadiationLengthInM();
	const double b_N_ref = tcol->GetMaterial()->GetSixtrackNuclearSlope();

	//We have now read the material properties, now to scale these if required to the current energy scale etc
	center_of_mass_squared = 2 * ProtonMassMeV * MeV * E0;  //ecmsq in SixTrack
	//pp cross-sections and parameters for energy dependence scaling
	const double p_reference = 450.0 * GeV;         //Reference energy at which the scattering data is based on.	(pref)

	const double pp_total_reference = 0.04;         //Proton total cross section at reference energy		(pptref)

	const double pp_elastic_reference = 0.007;      //Proton elastic cross section at reference energy		(pperef)

	const double single_diffractive_constant = 0.00068; //Single Diffractive constant					(sdcoe)

	const double pp_total_constant = 0.05788;       //pp total cross section constant				(pptco)

	const double pp_elastic_constant = 0.04792;     //pp elastic scattering cross section power constant		(ppeco)

	const double free_nucleon_constant = 1.618;     //free nucleon constant						(freeco)

	t_low_cut = 0.0009982;                  //Rutherford scattering cut scale (GeV^2)			(tlcut)
	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME
	//ppsd = sdcoe * log(0.15d0 * ecmsq)

	//Atom radius
	double atomic_radius = 1.2e-15 * pow(A, (1.0 / 3.0));  // In m, remember elsewhere if using area, to convert to barns

	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A, 1.0 / 3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.
	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy

	double sigma_pp_total = pp_total_reference * pow((E0 / p_reference), pp_total_constant);     //(pptot)

	//Next we also scale the elastic proton-nucleon cross section to the current energy
	//Remember this is for a single nucleon, when colplying to a nucleus
	const double sigma_pp_elastic = pp_elastic_reference * pow((E0 / p_reference), pp_elastic_constant); // demolaize equation 1.21 / catalan 3.7

	//And here the elastic cross section is adjusted to the number of nucleons in this material
	sigma_pn_elastic = free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	const double sigma_pp_SingleDiffractive = single_diffractive_constant * log(0.15 * center_of_mass_squared);

	//And again scale to the number of nucleons
	sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);

	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference * sigma_pN_total / sigma_pN_total_reference;

	//Caluclate the full nucleus elastic contribution
	sigma_pN_elastic = sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive;

	//Work on slopes next
	b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)); // slope given on GeV units

	b_N = b_N_ref * (sigma_pN_total / sigma_pN_total_reference);

	//Finally calculate the mean free path
	lambda_tot = A * 1.e-6 / ((sigma_pN_total + sigma_Rutherford) * barn * rho * Avogadro); // total mean free path (units meter)

	SetScatterConfigured(true);
	bool output_scattering_details = false;
	if(output_scattering_details)
	{
		std::cout << "SixTrack Scattering Configuration" << std::endl;
		std::cout << "..............................." << std::endl;
		std::cout << "Material : " << name << std::endl;
		std::cout << "Look at the material database reference:" << std::endl;
		std::cout << "Material reference proton-Nucleus total scattering cross section = "
				  << sigma_pN_total_reference << std::endl;
		std::cout << "Material reference proton-Nucleus inelastic scattering cross section = "
				  << sigma_pN_inelastic_reference << std::endl;
		std::cout << "Material reference Rutherford scattering cross section = " << sigma_Rutherford_reference
				  << std::endl;
		std::cout << "dEdx = " << dEdx << std::endl;
		std::cout << " rho = " << rho << std::endl;                     //density (g /cm^3)
		std::cout << "A = " << A << std::endl;
		std::cout << "Z = " << Z << std::endl;
		std::cout << "X0 = " << X0 << std::endl;
		std::cout << "b_N_ref = " << b_N_ref << std::endl;
		std::cout << "center_of_mass_squared = " << center_of_mass_squared << std::endl;
		std::cout << "Reference energy at which the scattering data is based on (pref) = " << p_reference << std::endl;
		std::cout << "Proton total cross section at reference energy (pptref) = " << pp_total_reference << std::endl;
		std::cout << "Proton elastic cross section at reference energy = " << pp_elastic_reference << std::endl;
		std::cout << "Single Diffractive constant = " << single_diffractive_constant << std::endl;
		std::cout << "pp total cross section constant = " << pp_total_constant << std::endl;
		std::cout << "pp elastic scattering cross section power constant = " << pp_elastic_constant << std::endl;
		std::cout << "Rutherford scattering cut scale (GeV^2) = " << t_low_cut << std::endl;
		std::cout << "atomic_radius = " << atomic_radius << std::endl;
		std::cout << "free nucleon constant = " << free_nucleon_constant << std::endl;
		std::cout << " number of free nucleons available for scattering= " << free_nucleon_count << std::endl;
		std::cout << "Scaled Elastic Cross section = " << sigma_pp_elastic << std::endl;
		std::cout << "Scaled tot cross section = " << sigma_pp_total << std::endl;
		std::cout << "Scaled Elastic Cross section Nucleus tot = " << sigma_pn_elastic << std::endl;
		std::cout << "Scaled Single diffractive cross section = " << sigma_pp_SingleDiffractive << std::endl;
		std::cout << "Scaled Single diffractive cross section Nucleus tot = " << sigma_pn_SingleDiffractive
				  << std::endl;
		std::cout << "sigma Rutherford(bug : not scaled)" << sigma_Rutherford << std::endl;
		std::cout << "Scaled total without Rutherford inclusion = " << sigma_pN_total << std::endl;
		std::cout << "Scaled Total inelastic = " << sigma_pN_inelastic << std::endl;
		std::cout << "Scaled full nucleus elastic contribution = " << sigma_pN_elastic << std::endl;
		std::cout << "b_pp(GeV units) = " << b_pp << std::endl;
		std::cout << "b_N (GeV units)= " << b_N << std::endl;
		std::cout << "Total mean free path (units meter) = " << lambda_tot << std::endl;
	}

}

int ProtonBunch::ScatterSixtrack(PSvector& p, double x, const Collimator* col)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference energy
	// col is the colerture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		OK
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

	if(!ScatterConfigured)
	{
		ConfigureScatter(col);
	}
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);

	//Keep track of distance along the collimator for colerture checking (colerture could vary with z)
	double z = 0;   // distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05; // maximum allowed energy loss - 95%

	while(x > 0)
	{
		bool interacted;    //Interacted on this step or not - true/false?
		double t = 0.0;       //Momentum transfer
		double delta_s = -lambda_tot * log(RandomNG::uniform(0, 1));
		double step_size;

		interacted = (x > delta_s);
		step_size = interacted ? delta_s : x;

		//Do MCS + dE/dx
		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());  // particle energy
		double thick = step_size / X0;  // material length in radiation lengths

		//Start dEdx
		double dp = dEdx * step_size;
		//End dEdx

		double E2 = E1 - dp;

		if(E2 <= 1.0)
		{
			p.ct() = z;
			return 1;
		}

		p.dp() = ((E1 - dp) - E0) / E0;
		double Eav = (E1 + E2) / 2.0;

		/*
		 *
		 *	MULTIPLE COULOMB SCATTERING
		 *
		 */
		double theta0 = 13.6 * MeV * sqrt(thick) * (1.0 + 0.038 * log(thick)) / Eav;    // small-angle Coulomb scattering

		pair<double, double> s = CoulombScatter(step_size, theta0);
		p.x() += s.first;
		p.xp() += s.second;

		s = CoulombScatter(step_size, theta0);
		p.y() += s.first;
		p.yp() += s.second;

		if(E2 < (E0 / 100.0))
		{
			return 4;
		}

		//Check we are still inside the collimator
		//If not, the particle leaves
		z += zstep;
		if(tcol->GetAperture()->CheckWithinApertureBoundaries(p.x(), p.y(), int_s + z))
		{
			tally[0]++;
			p.x() += p.xp() * x;
			p.y() += p.yp() * x;
			return 0;
		}

		//Reset this, not all interactions use it...
		dp = 0.0;

		//Point process interaction
		if(interacted)
		{
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0, 1) * (sigma_pN_total + sigma_Rutherford);

			//Choose which scattering process to do
			/*
			 *
			 *	Elastic scatter pN (proton - Nucleus)
			 *
			 */
			if((r -= sigma_pN_elastic) < 0)
			{
				tally[1]++;
				t = -log(RandomNG::uniform(0, 1)) / b_N;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Elastic scatter pn (proton - nucleon)
			 *
			 */
			else if((r -= sigma_pn_elastic) < 0)
			{
				tally[2]++;
				t = -log(RandomNG::uniform(0, 1)) / b_pp;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Single Diffractive
			 *
			 */

			else if((r -= sigma_pn_SingleDiffractive) < 0)
			{
				tally[3]++;

				double xm2 = exp(RandomNG::uniform(0, 1) * log(0.15 * center_of_mass_squared));
				double b = 0.0;
				if(xm2 < 2.0)
				{
					b = 2 * b_pp;
				}
				else if(2.0 <= xm2 && xm2 <= 5.0)
				{
					b = (106.0 - 17.0 * xm2) * b_pp / 26.0;    //This isn't what is in the sixtrack code (/26.0 typo (should be /36.0)) but is what is listed in N. Lasheras' thesis...
				}
				else if(xm2 > 5.0)
				{
					b = 7.0 * b_pp / 12.0;
				}
				t = -log(RandomNG::uniform(0, 1)) / b;
				dp = xm2 * E1 / center_of_mass_squared;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 2;
				p.sd() = 1;
			}

			/*
			 *
			 *	Rutherford coulomb scattering
			 *
			 */
			else if((r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut = t_low_cut;
				t = tcut / (1 - RandomNG::uniform(0, 1)); // generates 1/t squared distribution,
				p.type() = 3;
			}

			/*
			 *
			 *	Inelastic interaction - no more proton :(
			 *
			 */
			else
			{
				tally[5]++;
				p.ct() = z;
				return 1;
			}
			double E3 = E2 - dp;

			if(E3 <= 0.10)
			{
				p.ct() = z;
				return 1;
			}

			p.dp() = (E3 - E0) / E0;

			double theta = sqrt(t) / E3;
			double phi = RandomNG::uniform(-pi, pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
			if(std::isnan(p.xp()) || std::isnan(p.yp()))
			{
				cout << p << endl;
				cout << t << "\t" << theta << "\t" << E3 << endl;
			}
		}

		x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP || p.dp() < -1.0)    //dp cut should be at 95%
	{
		p.ct() = z;
		return 1;
	}

	return 0;
}

/****************************************************************
 *
 *
 *	"Sixtrack/K2 style with advanced ionization process" scattering physics
 *
 *
 ****************************************************************/

void ProtonBunch::ConfigureScatterSixtrackAdvancedIonization(const Collimator* col)
{
	//Do a cast to check if we have a "Collimator"
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);
	if(!tcol)
	{
		throw MerlinException("ScatterProton : No Collimator colerture");
	}

	double P0 = GetReferenceMomentum();
	E0 = sqrt(P0 * P0 + pow(ProtonMassMeV * MeV, 2));

	/**
	 * The proton gamma factor
	 */
	double gamma = E0 / (ProtonMassMeV * MeV);

	/**
	 * The proton beta
	 */
	double beta = sqrt(1 - (1 / (gamma * gamma)));

	const double sigma_pN_total_reference = tcol->GetMaterial()->GetSixtrackTotalNucleusCrossSection();      //Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tcol->GetMaterial()->GetSixtrackInelasticNucleusCrossSection();  //Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tcol->GetMaterial()->GetSixtrackRutherfordCrossSection();      //Material reference Rutherford scattering cross section
	rho = tcol->GetMaterial()->GetDensity() / 1000.0;
	A = tcol->GetMaterial()->GetAtomicMass();                        //Atomic mass
	Z = tcol->GetMaterial()->GetAtomicNumber();              //Atomic number
	X0 = tcol->GetMaterial()->GetRadiationLengthInM();
	I = tcol->GetMaterial()->GetMeanExcitationEnergy() / eV;
	const double ElectronDensity = tcol->GetMaterial()->GetElectronDensity();        //N_e / m^3
	const double PlasmaEnergy = tcol->GetMaterial()->GetPlasmaEnergy() / eV;
	const double b_N_ref = tcol->GetMaterial()->GetSixtrackNuclearSlope();
	//We have now read the material properties, now to scale these if required to the current energy scale etc
	center_of_mass_squared = 2 * ProtonMassMeV * MeV * E0;  //ecmsq in SixTrack
	//pp cross-sections and parameters for energy dependence scaling
	const double p_reference = 450.0 * GeV;         //Reference energy at which the scattering data is based on.	(pref)

	const double pp_total_reference = 0.04;         //Proton total cross section at reference energy		(pptref)

	const double pp_elastic_reference = 0.007;      //Proton elastic cross section at reference energy		(pperef)

	const double single_diffractive_constant = 0.00068; //Single Diffractive constant					(sdcoe)

	const double pp_total_constant = 0.05788;       //pp total cross section constant				(pptco)

	const double pp_elastic_constant = 0.04792;     //pp elastic scattering cross section power constant		(ppeco)

	const double free_nucleon_constant = 1.618;     //free nucleon constant						(freeco)

	t_low_cut = 0.0009982;                  //Rutherford scattering cut scale (GeV^2)			(tlcut)
	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME

	//Atom radius
	double atomic_radius = 1.2e-15 * pow(A, (1.0 / 3.0));  // In m, remember elsewhere if using area, to convert to barns

	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A, 1.0 / 3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.
	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy

	double sigma_pp_total = pp_total_reference * pow((E0 / p_reference), pp_total_constant);     //(pptot)

	//Next we also scale the elastic proton-nucleon cross section to the current energy
	//Remember this is for a single nucleon, when colplying to a nucleus
	const double sigma_pp_elastic = pp_elastic_reference * pow((E0 / p_reference), pp_elastic_constant); // demolaize equation 1.21 / catalan 3.7

	//And here the elastic cross section is adjusted to the number of nucleons in this material
	sigma_pn_elastic = free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	const double sigma_pp_SingleDiffractive = single_diffractive_constant * log(0.15 * center_of_mass_squared);

	//And again scale to the number of nucleons
	sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);

	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference * sigma_pN_total / sigma_pN_total_reference;

	//Caluclate the full nucleus elastic contribution
	sigma_pN_elastic = sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive;

	//Work on slopes next
	b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)); // slope given on GeV units

	b_N = b_N_ref * (sigma_pN_total / sigma_pN_total_reference);

	tmax = (2 * ElectronMassMeV * beta * beta * gamma * gamma) / (1 + (2 * gamma * (ElectronMassMeV / ProtonMassMeV))
		+ pow((ElectronMassMeV / ProtonMassMeV), 2)) * MeV;

	static const double xi1 = 2.0 * pi * pow(ElectronRadius, 2) * ElectronMass * pow(SpeedOfLight, 2);
	xi0 = xi1 * ElectronDensity;

	C = 1 + 2 * log(I / PlasmaEnergy);

	if((I / eV) < 100)
	{
		if(C <= 3.681)
		{
			C0 = 0.2;
			C1 = 2.0;
		}
		else
		{
			C0 = 0.326 * C - 1.0;
			C1 = 2.0;
		}
	}
	else    //I >= 100eV
	{
		if(C <= 5.215)
		{
			C0 = 0.2;
			C1 = 3.0;
		}
		else
		{
			C0 = 0.326 * C - 1.5;
			C1 = 3.0;
		}
	}

	//Finally calculate the mean free path
	lambda_tot = A * 1.e-6 / ((sigma_pN_total + sigma_Rutherford) * barn * rho * Avogadro); // total mean free path (units meter)
	SetScatterConfigured(true);

	bool output_scattering_details = false;
	if(output_scattering_details)
	{
		std::cout << "SixTrack with advanced ionization process configuration" << std::endl;
		std::cout << "Material reference proton-Nucleus total scattering cross section = "
				  << sigma_pN_total_reference << std::endl;
		std::cout << "Material reference proton-Nucleus inelastic scattering cross section = "
				  << sigma_pN_inelastic_reference << std::endl;
		std::cout << "Material reference Rutherford scattering cross section = " << sigma_Rutherford_reference
				  << std::endl;
		std::cout << "dEdx = " << dEdx << std::endl;
		std::cout << " rho = " << rho << std::endl;                     //density (g /cm^3)
		std::cout << "A = " << A << std::endl;
		std::cout << "Z = " << Z << std::endl;
		std::cout << "X0 = " << X0 << std::endl;
		std::cout << "b_N_ref = " << b_N_ref << std::endl;
		std::cout << "center_of_mass_squared = " << center_of_mass_squared << std::endl;
		std::cout << "Proton total cross section at reference energy (pptref) = " << pp_total_reference << std::endl;
		std::cout << "Reference energy at which the scattering data is based on (pref) = " << p_reference << std::endl;
		std::cout << "Single Diffractive constant = " << single_diffractive_constant << std::endl;
		std::cout << "Proton elastic cross section at reference energy = " << pp_elastic_reference << std::endl;
		std::cout << "pp total cross section constant = " << pp_total_constant << std::endl;
		std::cout << "pp elastic scattering cross section power constant = " << pp_elastic_constant << std::endl;
		std::cout << "free nucleon constant = " << free_nucleon_constant << std::endl;
		std::cout << "Rutherford scattering cut scale (GeV^2) = " << t_low_cut << std::endl;
		std::cout << "atomic_radius = " << atomic_radius << std::endl;
		std::cout << " number of free nucleons available for scattering= " << free_nucleon_count << std::endl;
		std::cout << "Scaled tot cross section = " << sigma_pp_total << std::endl;
		std::cout << "Scaled Elastic Cross section = " << sigma_pp_elastic << std::endl;
		std::cout << "Scaled Elastic Cross section Nucleus tot = " << sigma_pn_elastic << std::endl;
		std::cout << "Scaled Single diffractive cross section = " << sigma_pp_SingleDiffractive << std::endl;
		std::cout << "Scaled Single diffractive cross section Nucleus tot = " << sigma_pn_SingleDiffractive
				  << std::endl;
		std::cout << "sigma Rutherford(bug : not scaled)" << sigma_Rutherford << std::endl;
		std::cout << "Scaled total without Rutherford inclusion = " << sigma_pN_total << std::endl;
		std::cout << "Scaled Total inelastic = " << sigma_pN_inelastic << std::endl;
		std::cout << "Scaled full nucleus elastic contribution = " << sigma_pN_elastic << std::endl;
		std::cout << "b_pp(GeV units) = " << b_pp << std::endl;
		std::cout << "b_N (GeV units)= " << b_N << std::endl;
		std::cout << "Tmax: " << tmax / GeV << " GeV" << std::endl;
		std::cout << "Total mean free path (units meter) = " << lambda_tot << std::endl;
	}

}

int ProtonBunch::ScatterSixtrackAdvancedIonization(PSvector& p, double x, const Collimator* col)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference energy
	// col is the colerture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		OK
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

	if(!ScatterConfigured)
	{
		ConfigureScatter(col);
	}
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);

	//Keep track of distance along the collimator for colerture checking (colerture could vary with z)
	double z = 0;   // distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05; // maximum allowed energy loss - 95%

	while(x > 0)
	{
		bool interacted;    //Interacted on this step or not - true/false?
		double t = 0.0;       //Momentum transfer
		double delta_s = -lambda_tot * log(RandomNG::uniform(0, 1));
		double step_size;

		interacted = (x > delta_s);
		step_size = interacted ? delta_s : x;

		//Do MCS + dE/dx
		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());  // particle energy
		double thick = step_size / X0;  // material length in radiation lengths

		//Start dEdx
		double gamma = E1 / (ProtonMassMeV * MeV);
		double beta = sqrt(1 - (1 / (gamma * gamma)));
		double land = RandomNG::landau();

		double xi = (xi0 * step_size / (beta * beta)) / ElectronCharge * (eV / MeV);

		//Density correction
		double ddx = log10(beta * gamma);
		if(ddx > C1)
		{
			delta = 4.606 * ddx - C;
		}
		else if(ddx >= C0 && ddx <= C1)
		{
			double m = 3.0;
			double xa = C / 4.606;
			double a = 4.606 * (xa - C0) / pow((C1 - C0), m);
			delta = 4.606 * ddx - C + a * pow((C1 - ddx), m);
		}
		else
		{
			delta = 0.0;
		}

		double tcut = 2.0 * MeV;
		tcut = tmax;

		//Mott Correction
		double G = pi * FineStructureConstant * beta / 2.0;
		double q = (2 * (tmax / MeV) * (ElectronMassMeV)) / (pow((0.843 / MeV), 2));
		double S = log(1 + q);

		double L1 = 0.0;

		double yL2 = FineStructureConstant / beta;

		double L2sum = 1.202001688211;  //Sequence limit calculated with mathematica
		double L2 = -yL2 * yL2 * L2sum;

		double F = G - S + 2 * (L1 + L2);
		double deltaE = xi * (log(2 * ElectronMassMeV * beta * beta * gamma * gamma * (tcut / MeV) / pow(I / MeV, 2))
			- (beta * beta) * (1 + ((tcut / MeV) / (tmax / MeV))) - delta + F - 1.0 - euler);
		deltaE = xi * (log(2 * ElectronMassMeV * beta * beta * gamma * gamma * xi / pow(I / MeV, 2)) - (beta * beta)
			- delta + F + 0.20);
		if(xi > (2 * ElectronMassMeV * beta * beta * gamma * gamma))
		{
			cout << "de,xi " << deltaE << "\t" << xi << endl;
		}

		double dp = ((xi * land) - deltaE) * MeV;
		//End dEdx

		double E2 = E1 - dp;

		if(E2 <= 1.0)
		{
			p.ct() = z;
			return 1;
		}

		p.dp() = ((E1 - dp) - E0) / E0;
		double Eav = (E1 + E2) / 2.0;

		/*
		 *
		 *	MULTIPLE COULOMB SCATTERING
		 *
		 */
		double theta0 = 13.6 * MeV * sqrt(thick) * (1.0 + 0.038 * log(thick)) / Eav;    // small-angle Coulomb scattering

		pair<double, double> s = CoulombScatter(step_size, theta0);
		p.x() += s.first;
		p.xp() += s.second;

		s = CoulombScatter(step_size, theta0);
		p.y() += s.first;
		p.yp() += s.second;

		if(E2 < (E0 / 100.0))
		{
			return 4;
		}

		//Check we are still inside the collimator
		//If not, the particle leaves
		z += zstep;
		if(tcol->GetAperture()->CheckWithinApertureBoundaries(p.x(), p.y(), int_s + z))
		{
			tally[0]++;
			p.x() += p.xp() * x;
			p.y() += p.yp() * x;
			return 0;
		}

		//Reset this, not all interactions use it...
		dp = 0.0;

		//Point process interaction
		if(interacted)
		{
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0, 1) * (sigma_pN_total + sigma_Rutherford);

			//Choose which scattering process to do
			/*
			 *
			 *	Elastic scatter pN (proton - Nucleus)
			 *
			 */
			if((r -= sigma_pN_elastic) < 0)
			{
				tally[1]++;
				t = -log(RandomNG::uniform(0, 1)) / b_N;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Elastic scatter pn (proton - nucleon)
			 *
			 */
			else if((r -= sigma_pn_elastic) < 0)
			{
				tally[2]++;
				t = -log(RandomNG::uniform(0, 1)) / b_pp;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Single Diffractive
			 *
			 */

			else if((r -= sigma_pn_SingleDiffractive) < 0)
			{
				tally[3]++;

				double xm2 = exp(RandomNG::uniform(0, 1) * log(0.15 * center_of_mass_squared));
				double b = 0.0;
				if(xm2 < 2.0)
				{
					b = 2 * b_pp;
				}
				else if(2.0 <= xm2 && xm2 <= 5.0)
				{
					b = (106.0 - 17.0 * xm2) * b_pp / 26.0;    //This isn't what is in the sixtrack code (/26.0 typo (should be /36.0)) but is what is listed in N. Lasheras' thesis...
				}
				else if(xm2 > 5.0)
				{
					b = 7.0 * b_pp / 12.0;
				}
				t = -log(RandomNG::uniform(0, 1)) / b;
				dp = xm2 * E1 / center_of_mass_squared;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;

				p.type() = 2;
				p.sd() = 1;
			}

			/*
			 *
			 *	Rutherford coulomb scattering
			 *
			 */
			else if((r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut = t_low_cut;
				t = tcut / (1 - RandomNG::uniform(0, 1)); // generates 1/t squared distribution,
				p.type() = 3;
			}

			/*
			 *
			 *	Inelastic interaction - no more proton :(
			 *
			 */
			else
			{
				tally[5]++;
				p.ct() = z;
				return 1;
			}
			double E3 = E2 - dp;

			if(E3 <= 0.10)
			{
				p.ct() = z;
				return 1;
			}

			p.dp() = (E3 - E0) / E0;

			double theta = sqrt(t) / E3;
			double phi = RandomNG::uniform(-pi, pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
			if(std::isnan(p.xp()) || std::isnan(p.yp()))
			{
				cout << p << endl;
				cout << t << "\t" << theta << "\t" << E3 << endl;
			}
		}

		x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP || p.dp() < -1.0)    //dp cut should be at 95%
	{
		p.ct() = z;
		return 1;
	}

	return 0;
}

/****************************************************************
 *
 *
 *	"Sixtrack/K2 style with advanced pomeron elastic scattering" scattering physics
 *
 *
 ****************************************************************/

void ProtonBunch::ConfigureScatterSixtrackAdvancedElastic(const Collimator* col)
{
	//Do a cast to check if we have a "Collimator"
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);
	if(!tcol)
	{
		throw MerlinException("ScatterProton : No Collimator colerture");
	}

	double P0 = GetReferenceMomentum();
	E0 = sqrt(P0 * P0 + pow(ProtonMassMeV * MeV, 2));

	if(!GotElastic)
	{
		ElasticScatter = new ppElasticScatter();
		/**
		 * Generate the elastic differential cross section for the proton energy
		 */
		ElasticScatter->SetTMin(1e-4);
		ElasticScatter->SetTMax(1.0);
		ElasticScatter->SetStepSize(1e-4);
		ElasticScatter->GenerateTDistribution(E0);
		GotElastic = true;
	}

	const double sigma_pN_total_reference = tcol->GetMaterial()->GetSixtrackTotalNucleusCrossSection();      //Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tcol->GetMaterial()->GetSixtrackInelasticNucleusCrossSection();  //Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tcol->GetMaterial()->GetSixtrackRutherfordCrossSection();      //Material reference Rutherford scattering cross section
	dEdx = tcol->GetMaterial()->GetSixtrackdEdx();                   //dE/dx - (GeV m^-1)
	rho = tcol->GetMaterial()->GetDensity() / 1000.0;
	A = tcol->GetMaterial()->GetAtomicMass();                        //Atomic mass
	Z = tcol->GetMaterial()->GetAtomicNumber();              //Atomic number
	X0 = tcol->GetMaterial()->GetRadiationLengthInM();

	const double b_N_ref = tcol->GetMaterial()->GetSixtrackNuclearSlope();

	//We have now read the material properties, now to scale these if required to the current energy scale etc
	center_of_mass_squared = 2 * ProtonMassMeV * MeV * E0;  //ecmsq in SixTrack
	//pp cross-sections and parameters for energy dependence scaling
	const double p_reference = 450.0 * GeV;         //Reference energy at which the scattering data is based on.	(pref)

	const double pp_total_reference = 0.04;         //Proton total cross section at reference energy		(pptref)

	const double single_diffractive_constant = 0.00068; //Single Diffractive constant					(sdcoe)

	const double pp_total_constant = 0.05788;       //pp total cross section constant				(pptco)

	const double free_nucleon_constant = 1.618;     //free nucleon constant						(freeco)

	t_low_cut = 0.0009982;                  //Rutherford scattering cut scale (GeV^2)			(tlcut)

	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME
	//ppsd = sdcoe * log(0.15d0 * ecmsq)

	//Atom radius
	double atomic_radius = 1.2e-15 * pow(A, (1.0 / 3.0));  // In m, remember elsewhere if using area, to convert to barns

	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A, 1.0 / 3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.
	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy

	double sigma_pp_total = pp_total_reference * pow((E0 / p_reference), pp_total_constant);     //(pptot)

	/**
	 * The total elastic cross section
	 */
	const double sigma_pp_elastic = ElasticScatter->GetElasticCrossSectionN();
	const double sigma_pp_elasticEM = ElasticScatter->GetElasticCrossSection();
	double ElasticDifference = sigma_pp_elasticEM - sigma_pp_elastic;

	//And here the elastic cross section is adjusted to the number of nucleons in this material
	sigma_pn_elastic = free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	const double sigma_pp_SingleDiffractive = single_diffractive_constant * log(0.15 * center_of_mass_squared);

	//And again scale to the number of nucleons
	sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);

	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference * sigma_pN_total / sigma_pN_total_reference;

	//Caluclate the full nucleus elastic contribution
	sigma_pN_elastic = sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive;

	//Work on slopes next
	b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)); // slope given on GeV units

	b_N = b_N_ref * (sigma_pN_total / sigma_pN_total_reference);

	//Finally calculate the mean free path
	lambda_tot = A * 1.e-6 / ((sigma_pN_total + sigma_Rutherford + ElasticDifference) * barn * rho * Avogadro); // total mean free path (units meter)

	SetScatterConfigured(true);
	bool output_scattering_details = false;
	if(output_scattering_details)
	{
		std::cout << "SixTrack with advanced Elastic Configuration" << std::endl;
		std::cout << "Material reference proton-Nucleus total scattering cross section = "
				  << sigma_pN_total_reference << std::endl;
		std::cout << "Material reference proton-Nucleus inelastic scattering cross section = "
				  << sigma_pN_inelastic_reference << std::endl;
		std::cout << "Material reference Rutherford scattering cross section = " << sigma_Rutherford_reference
				  << std::endl;
		std::cout << "dEdx = " << dEdx << std::endl;
		std::cout << " rho = " << rho << std::endl;                     //density (g /cm^3)
		std::cout << "A = " << A << std::endl;
		std::cout << "Z = " << Z << std::endl;
		std::cout << "X0 = " << X0 << std::endl;
		std::cout << "b_N_ref = " << b_N_ref << std::endl;
		std::cout << "center_of_mass_squared = " << center_of_mass_squared << std::endl;
		std::cout << "Reference energy at which the scattering data is based on (pref) = " << p_reference << std::endl;
		std::cout << "Proton total cross section at reference energy (pptref) = " << pp_total_reference << std::endl;
		std::cout << "Single Diffractive constant = " << single_diffractive_constant << std::endl;
		std::cout << "pp total cross section constant = " << pp_total_constant << std::endl;
		std::cout << "free nucleon constant = " << free_nucleon_constant << std::endl;
		std::cout << "Rutherford scattering cut scale (GeV^2) = " << t_low_cut << std::endl;
		std::cout << "atomic_radius = " << atomic_radius << std::endl;
		std::cout << " number of free nucleons available for scattering= " << free_nucleon_count << std::endl;
		std::cout << "Scaled tot cross section = " << sigma_pp_total << std::endl;
		std::cout << "Elastic pp cross section = " << sigma_pp_elastic << std::endl;
		std::cout << "Elastic pp cross section EM = " << sigma_pp_elasticEM << std::endl;
		std::cout << "Elastic difference = " << ElasticDifference << std::endl;
		std::cout << "Scaled Elastic Cross section nucleons tot = " << sigma_pn_elastic << std::endl;
		std::cout << "Scaled Single diffractive cross section = " << sigma_pp_SingleDiffractive << std::endl;
		std::cout << "Scaled Single diffractive cross section Nucleus tot = " << sigma_pn_SingleDiffractive
				  << std::endl;
		std::cout << "sigma Rutherford(bug : not scaled)" << sigma_Rutherford << std::endl;
		std::cout << "Scaled total without Rutherford inclusion = " << sigma_pN_total << std::endl;
		std::cout << "Scaled Total inelastic = " << sigma_pN_inelastic << std::endl;
		std::cout << "Scaled full nucleus elastic contribution = " << sigma_pN_elastic << std::endl;
		std::cout << "b_pp(GeV units) = " << b_pp << std::endl;
		std::cout << "b_N (GeV units)= " << b_N << std::endl;
		std::cout << "Total mean free path (units meter) = " << lambda_tot << std::endl;
	}

}

int ProtonBunch::ScatterSixtrackAdvancedElastic(PSvector& p, double x, const Collimator* col)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference energy
	// col is the colerture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		OK
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

	if(!ScatterConfigured)
	{
		ConfigureScatter(col);
	}
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);

	//Keep track of distance along the collimator for colerture checking (colerture could vary with z)
	double z = 0;   // distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05; // maximum allowed energy loss - 95%

	while(x > 0)
	{
		bool interacted;    //Interacted on this step or not - true/false?
		double t = 0.0;       //Momentum transfer
		double delta_s = -lambda_tot * log(RandomNG::uniform(0, 1));
		double step_size;

		interacted = (x > delta_s);
		step_size = interacted ? delta_s : x;

		//Do MCS + dE/dx
		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());  // particle energy
		double thick = step_size / X0;  // material length in radiation lengths

		//Start dEdx
		double dp = dEdx * step_size;
		//End dEdx

		double E2 = E1 - dp;

		if(E2 <= 1.0)
		{
			p.ct() = z;
			return 1;
		}

		p.dp() = ((E1 - dp) - E0) / E0;
		double Eav = (E1 + E2) / 2.0;

		/*
		 *
		 *	MULTIPLE COULOMB SCATTERING
		 *
		 */
		double theta0 = 13.6 * MeV * sqrt(thick) * (1.0 + 0.038 * log(thick)) / Eav;    // small-angle Coulomb scattering

		pair<double, double> s = CoulombScatter(step_size, theta0);
		p.x() += s.first;
		p.xp() += s.second;

		s = CoulombScatter(step_size, theta0);
		p.y() += s.first;
		p.yp() += s.second;

		if(E2 < (E0 / 100.0))
		{
			return 4;
		}

		//Check we are still inside the collimator
		//If not, the particle leaves
		z += zstep;
		if(tcol->GetAperture()->CheckWithinApertureBoundaries(p.x(), p.y(), int_s + z))
		{
			tally[0]++;
			p.x() += p.xp() * x;
			p.y() += p.yp() * x;
			return 0;
		}

		//Reset this, not all interactions use it...
		dp = 0.0;

		//Point process interaction
		if(interacted)
		{
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0, 1) * (sigma_pN_total + sigma_Rutherford);

			//Choose which scattering process to do
			/*
			 *
			 *	Elastic scatter pN (proton - Nucleus)
			 *
			 */
			if((r -= sigma_pN_elastic) < 0)
			{
				tally[1]++;
				t = -log(RandomNG::uniform(0, 1)) / b_N;
				double TargetMass = A * AtomicMassUnit; // Nucleus mass
				dp = t / (2 * TargetMass);
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Elastic scatter pn (proton - nucleon)
			 *
			 */
			else if((r -= sigma_pn_elastic) < 0)
			{
				tally[2]++;
				t = ElasticScatter->SelectT();
				dp = t / (2 * AtomicMassUnit);
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Single Diffractive
			 *
			 */

			else if((r -= sigma_pn_SingleDiffractive) < 0)
			{
				tally[3]++;

				double xm2 = exp(RandomNG::uniform(0, 1) * log(0.15 * center_of_mass_squared));
				double b = 0.0;
				if(xm2 < 2.0)
				{
					b = 2 * b_pp;
				}
				else if(2.0 <= xm2 && xm2 <= 5.0)
				{
					b = (106.0 - 17.0 * xm2) * b_pp / 26.0;    //This isn't what is in the sixtrack code (/26.0 typo (should be /36.0)) but is what is listed in N. Lasheras' thesis...
				}
				else if(xm2 > 5.0)
				{
					b = 7.0 * b_pp / 12.0;
				}
				t = -log(RandomNG::uniform(0, 1)) / b;
				dp = xm2 * E1 / center_of_mass_squared;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 2;
				p.sd() = 1;
			}

			/*
			 *
			 *	Rutherford coulomb scattering
			 *
			 */
			else if((r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut = t_low_cut;
				t = tcut / (1 - RandomNG::uniform(0, 1)); // generates 1/t squared distribution,
				p.type() = 3;
			}

			/*
			 *
			 *	Inelastic interaction - no more proton :(
			 *
			 */
			else
			{
				tally[5]++;
				p.ct() = z;
				return 1;
			}
			double E3 = E2 - dp;

			if(E3 <= 0.10)
			{
				p.ct() = z;
				return 1;
			}

			p.dp() = (E3 - E0) / E0;

			double theta = sqrt(t) / E3;
			double phi = RandomNG::uniform(-pi, pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
			if(std::isnan(p.xp()) || std::isnan(p.yp()))
			{
				cout << p << endl;
				cout << t << "\t" << theta << "\t" << E3 << endl;
			}
		}

		x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP || p.dp() < -1.0)    //dp cut should be at 95%
	{
		p.ct() = z;
		return 1;
	}

	return 0;
}

/****************************************************************
 *
 *
 *	"Sixtrack/K2 style with advanced pomeron Single Diffraction scattering" scattering physics
 *
 *
 ****************************************************************/

void ProtonBunch::ConfigureScatterSixtrackAdvancedSingleDiffraction(const Collimator* col)
{
	//Do a cast to check if we have a "Collimator"
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);
	if(!tcol)
	{
		throw MerlinException("ScatterProton : No Collimator colerture");
	}

	double P0 = GetReferenceMomentum();
	E0 = sqrt(P0 * P0 + pow(ProtonMassMeV * MeV, 2));

	if(!GotDiffractive)
	{
		double sLHC = pow(114.6192348986088, 2);
		double Mproton = 0.938272013;
		double Mpion = 0.1349766;
		double Mmin2 = pow(Mproton + Mpion, 2);
		double xi_th = Mmin2 / sLHC; // (M_p + M_pion)^2/s

		DiffractiveScatter = new ppDiffractiveScatter();
		/**
		 * Generate the Single Diffractive differential cross section for the proton energy
		 */
		DiffractiveScatter->SetTMin(0.0001);
		DiffractiveScatter->SetTMax(4);
		DiffractiveScatter->SetTStepSize(1e-4);
		DiffractiveScatter->SetXiMin(xi_th); //Threshould at (M_proton + M_pion)^2/s
		DiffractiveScatter->SetXiMax(0.12);
		DiffractiveScatter->SetXiStepSize(1e-6);

		DiffractiveScatter->GenerateDistribution(E0);
		GotDiffractive = true;
	}

	const double sigma_pN_total_reference = tcol->GetMaterial()->GetSixtrackTotalNucleusCrossSection();      //Material reference proton-Nucleus total scattering cross section
	const double sigma_pN_inelastic_reference = tcol->GetMaterial()->GetSixtrackInelasticNucleusCrossSection();  //Material reference proton-Nucleus inelastic scattering cross section
	const double sigma_Rutherford_reference = tcol->GetMaterial()->GetSixtrackRutherfordCrossSection();      //Material reference Rutherford scattering cross section
	dEdx = tcol->GetMaterial()->GetSixtrackdEdx();                   //dE/dx - (GeV m^-1)
	rho = tcol->GetMaterial()->GetDensity() / 1000.0;
	A = tcol->GetMaterial()->GetAtomicMass();                        //Atomic mass
	Z = tcol->GetMaterial()->GetAtomicNumber();              //Atomic number
	X0 = tcol->GetMaterial()->GetRadiationLengthInM();
	const double b_N_ref = tcol->GetMaterial()->GetSixtrackNuclearSlope();

	//We have now read the material properties, now to scale these if required to the current energy scale etc
	center_of_mass_squared = 2 * ProtonMassMeV * MeV * E0;  //ecmsq in SixTrack
	//pp cross-sections and parameters for energy dependence scaling
	const double p_reference = 450.0 * GeV;         //Reference energy at which the scattering data is based on.	(pref)

	const double pp_total_reference = 0.04;         //Proton total cross section at reference energy		(pptref)

	const double pp_elastic_reference = 0.007;      //Proton elastic cross section at reference energy		(pperef)

	const double pp_total_constant = 0.05788;       //pp total cross section constant				(pptco)

	const double pp_elastic_constant = 0.04792;     //pp elastic scattering cross section power constant		(ppeco)

	const double free_nucleon_constant = 1.618;     //free nucleon constant						(freeco)

	t_low_cut = 0.0009982;                  //Rutherford scattering cut scale (GeV^2)			(tlcut)
	//emr - material atomic size for Rutherford - see NCL thesis for correct term FIXME

	//Atom radius
	double atomic_radius = 1.2e-15 * pow(A, (1.0 / 3.0));  // In m, remember elsewhere if using area, to convert to barns

	//First calculate the number of "free nucleons" available for scattering
	const double free_nucleon_count = free_nucleon_constant * pow(A, 1.0 / 3.0);

	//Could put the following block within the collimator iteration loop if one is being ultra-pedantic.
	//Cross sections need scaling from the reference energy to the beam energy
	//Nucleus cross sections need scaling to the number of nucleons
	//First task is to calculate the adjusted total cross section at this energy from the reference energy

	double sigma_pp_total = pp_total_reference * pow((E0 / p_reference), pp_total_constant);     //(pptot)

	//Next we also scale the elastic proton-nucleon cross section to the current energy
	//Remember this is for a single nucleon, when colplying to a nucleus
	const double sigma_pp_elastic = pp_elastic_reference * pow((E0 / p_reference), pp_elastic_constant); // demolaize equation 1.21 / catalan 3.7

	//And here the elastic cross section is adjusted to the number of nucleons in this material
	sigma_pn_elastic = free_nucleon_count * sigma_pp_elastic;

	//Single diffractive cross section
	double sigma_pp_SingleDiffractive = DiffractiveScatter->GetDiffractiveCrossSection();

	//And again scale to the number of nucleons
	sigma_pn_SingleDiffractive = free_nucleon_count * sigma_pp_SingleDiffractive;

	//Next fix Rutherford coulomb scattering
	//FIXME
	sigma_Rutherford = sigma_Rutherford_reference;

	//Correct total without Rutherford inclusion
	sigma_pN_total = sigma_pN_total_reference + free_nucleon_count * (sigma_pp_total - pp_total_reference);

	//And the inelastic
	double sigma_pN_inelastic = sigma_pN_inelastic_reference * sigma_pN_total / sigma_pN_total_reference;

	//Caluclate the full nucleus elastic contribution
	sigma_pN_elastic = sigma_pN_total - sigma_pN_inelastic - sigma_pn_elastic - sigma_pn_SingleDiffractive;

	//Work on slopes next
	b_pp = 8.5 + 1.086 * log(sqrt(center_of_mass_squared)); // slope given on GeV units

	b_N = b_N_ref * (sigma_pN_total / sigma_pN_total_reference);

	//Finally calculate the mean free path
	lambda_tot = A * 1.e-6 / ((sigma_pN_total + sigma_Rutherford) * barn * rho * Avogadro); // total mean free path (units meter)

	SetScatterConfigured(true);
	bool output_scattering_details = false;
	if(output_scattering_details)
	{
		std::cout << "Material reference proton-Nucleus total scattering cross section = "
				  << sigma_pN_total_reference << std::endl;
		std::cout << "Material reference proton-Nucleus inelastic scattering cross section = "
				  << sigma_pN_inelastic_reference << std::endl;
		std::cout << "Material reference Rutherford scattering cross section = " << sigma_Rutherford_reference
				  << std::endl;
		std::cout << "dEdx = " << dEdx << std::endl;
		std::cout << " rho = " << rho << std::endl;                     // density (g /cm^3)
		std::cout << "A = " << A << std::endl;
		std::cout << "Z = " << Z << std::endl;
		std::cout << "X0 = " << X0 << std::endl;
		std::cout << "b_N_ref = " << b_N_ref << std::endl;
		std::cout << "center_of_mass_squared = " << center_of_mass_squared << std::endl;
		std::cout << "Reference energy at which the scattering data is based on (pref) = " << p_reference << std::endl;
		std::cout << "Proton total cross section at reference energy (pptref) = " << pp_total_reference << std::endl;
		std::cout << "Proton elastic cross section at reference energy = " << pp_elastic_reference << std::endl;
		std::cout << "pp total cross section constant = " << pp_total_constant << std::endl;
		std::cout << "pp elastic scattering cross section power constant = " << pp_elastic_constant << std::endl;
		std::cout << "free nucleon constant = " << free_nucleon_constant << std::endl;
		std::cout << "Rutherford scattering cut scale (GeV^2) = " << t_low_cut << std::endl;
		std::cout << "atomic_radius = " << atomic_radius << std::endl;
		std::cout << " number of free nucleons available for scattering= " << free_nucleon_count << std::endl;
		std::cout << "Scaled tot cross section = " << sigma_pp_total << std::endl;
		std::cout << "Scaled Elastic Cross section = " << sigma_pp_elastic << std::endl;
		std::cout << "Scaled Elastic Cross section Nucleus tot = " << sigma_pn_elastic << std::endl;
		std::cout << "Scaled Single diffractive cross section Nucleus tot = " << sigma_pn_SingleDiffractive
				  << std::endl;
		std::cout << "sigma Rutherford(bug : not scaled)" << sigma_Rutherford << std::endl;
		std::cout << "Scaled total without Rutherford inclusion = " << sigma_pN_total << std::endl;
		std::cout << "Scaled Total inelastic = " << sigma_pN_inelastic << std::endl;
		std::cout << "Scaled full nucleus elastic contribution = " << sigma_pN_elastic << std::endl;
		std::cout << "b_pp(GeV units) = " << b_pp << std::endl;
		std::cout << "b_N (GeV units)= " << b_N << std::endl;
		std::cout << "Total mean free path (units meter) = " << lambda_tot << std::endl;
	}
}

int ProtonBunch::ScatterSixtrackAdvancedSingleDiffraction(PSvector& p, double x, const Collimator* col)
{
	// p is the scattering Proton - a single particle.
	// x is the length of the collimator
	// E0 is the reference energy
	// col is the colerture

	//We simulate many physical processes that occur in bulk materials:
	//Ionization	(-dE/dx)		TODO
	//multiple coulomb scattering		TODO

	//proton-nucleon elastic		OK
	//proton-nucleon inelastic		OK
	//proton-nucleus elastic		OK
	//proton-nucleus inelastic		OK

	//Single Diffractive			OK - CHECK 1/26 FACTOR
	//Rutherford scattering			TODO

	if(!ScatterConfigured)
	{
		ConfigureScatter(col);
	}
	const Collimator* tcol = dynamic_cast<const Collimator*> (col);

	//Keep track of distance along the collimator for colerture checking (colerture could vary with z)
	double z = 0;   // distance travelled along the collimator. Units (m)

	static const double MAXDP = 1.0 - 0.05; // maximum allowed energy loss - 95%

	while(x > 0)
	{
		bool interacted;    //Interacted on this step or not - true/false?
		double t = 0.0;       //Momentum transfer
		double delta_s = -lambda_tot * log(RandomNG::uniform(0, 1));
		double step_size;

		interacted = (x > delta_s);
		step_size = interacted ? delta_s : x;

		//Do MCS + dE/dx
		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		double E1 = E0 * (1 + p.dp());  // particle energy
		double thick = step_size / X0;  // material length in radiation lengths

		//Start dEdx
		double dp = dEdx * step_size;
		//End dEdx

		double E2 = E1 - dp;

		if(E2 <= 1.0)
		{
			p.ct() = z;
			return 1;
		}

		p.dp() = ((E1 - dp) - E0) / E0;
		double Eav = (E1 + E2) / 2.0;

		/*
		 *
		 *	MULTIPLE COULOMB SCATTERING
		 *
		 */
		double theta0 = 13.6 * MeV * sqrt(thick) * (1.0 + 0.038 * log(thick)) / Eav;    // small-angle Coulomb scattering

		pair<double, double> s = CoulombScatter(step_size, theta0);
		p.x() += s.first;
		p.xp() += s.second;

		s = CoulombScatter(step_size, theta0);
		p.y() += s.first;
		p.yp() += s.second;

		if(E2 < (E0 / 100.0))
		{
			return 4;
		}

		//Check we are still inside the collimator
		//If not, the particle leaves
		z += zstep;
		if(tcol->GetAperture()->CheckWithinApertureBoundaries(p.x(), p.y(), int_s + z))
		{
			tally[0]++;
			p.x() += p.xp() * x;
			p.y() += p.yp() * x;
			return 0;
		}

		//Reset this, not all interactions use it...
		dp = 0.0;

		//Point process interaction
		if(interacted)
		{
			E1 = E0 * (1 + p.dp());
			double r = RandomNG::uniform(0, 1) * (sigma_pN_total + sigma_Rutherford);

			//Choose which scattering process to do
			/*
			 *
			 *	Elastic scatter pN (proton - Nucleus)
			 *
			 */
			if((r -= sigma_pN_elastic) < 0)
			{
				tally[1]++;
				t = -log(RandomNG::uniform(0, 1)) / b_N;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Elastic scatter pn (proton - nucleon)
			 *
			 */
			else if((r -= sigma_pn_elastic) < 0)
			{
				tally[2]++;
				t = -log(RandomNG::uniform(0, 1)) / b_pp;
				double P0 = GetReferenceMomentum();
				double theta = sqrt(t) / P0;
				p.type() = 1;
			}

			/*
			 *
			 *	Single Diffractive
			 *
			 */

			else if((r -= sigma_pn_SingleDiffractive) < 0)
			{
				tally[3]++;
				std::pair<double, double> TM = DiffractiveScatter->Select();
				t = TM.first;
				double mrec = TM.second;

				dp = mrec * mrec * E1 / center_of_mass_squared;
				p.type() = 2;
				p.sd() = 1;
			}

			/*
			 *
			 *	Rutherford coulomb scattering
			 *
			 */
			else if((r -= sigma_Rutherford) < 0)
			{
				tally[4]++;
				double tcut = t_low_cut;
				t = tcut / (1 - RandomNG::uniform(0, 1)); // generates 1/t squared distribution,
				p.type() = 3;
			}

			/*
			 *
			 *	Inelastic interaction - no more proton :(
			 *
			 */
			else
			{
				tally[5]++;
				p.ct() = z;
				return 1;
			}
			double E3 = E2 - dp;

			if(E3 <= 0.10)
			{
				p.ct() = z;
				return 1;
			}

			p.dp() = (E3 - E0) / E0;

			double theta = acos(1 - (2 * pow(ProtonMassMeV * MeV, 2) - t) / (2 * E2 * E3));
			double phi = RandomNG::uniform(-pi, pi);
			p.xp() += theta * cos(phi);
			p.yp() += theta * sin(phi);
			if(std::isnan(p.xp()) || std::isnan(p.yp()))
			{
				cout << p << endl;
				cout << t << "\t" << theta << "\t" << E3 << endl;
			}
		}

		x -= step_size;
	} // end of while loop

	if(p.dp() < -MAXDP || p.dp() < -1.0)    //dp cut should be at 95%
	{
		p.ct() = z;
		return 1;
	}

	return 0;
}

void ProtonBunch::SetUpProfiling() const
{
	MERLIN_PROFILE_ADD_PROCESS("ProtonBunch::ConfigureScatterMerlin");
	MERLIN_PROFILE_ADD_PROCESS("ProtonBunch::ScatterMerlin");
}
