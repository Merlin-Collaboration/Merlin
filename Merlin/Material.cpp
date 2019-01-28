/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <iostream>

#include "Material.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

/*
 * Basic constructor taking the following arguments:
 * Name, Symbol, Atomic Mass, Atomic Number, Sigma_E, Sigma_I, Sigma_R, dEdx, Radiation Length, Density, Conductivity.
 */
Material::Material(std::string name0, std::string sym0, double A0, int AtomicNumber0, double Sigma_E0, double Sigma_I0,
	double Sigma_R0, double dEdx0, double X00, double Density0, double Conductivity0)
{
	Name = name0;
	Symbol = sym0;
	AtomicMass = A0;
	AtomicNumber = AtomicNumber0;
	sigma_pN_elastic = Sigma_E0;
	sigma_pN_inelastic = Sigma_I0;
	sigma_Rutherford = Sigma_R0;
	dEdx = dEdx0;
	X0 = X00;
	Density = Density0;
	Conductivity = Conductivity0;
}

double Material::CalculateElectronDensity()
{
	return AtomicNumber * Avogadro * Density * 0.001 / (AtomicMass * pow(centimeter, 3)); // n_e m^-3 (1e6 conversion from cm^3)
}

//returns the Plasma energy in GeV
double Material::CalculatePlasmaEnergy()
{
	return (PlanckConstantBar * sqrt((ElectronDensity * pow(ElectronCharge, 2)) / (ElectronMass
		   * FreeSpacePermittivity))) / ElectronCharge * eV;
}

double Material::CalculateMeanExcitationEnergy()
{
	return AtomicNumber * 10.0 * eV;
}

double Material::CalculateSixtrackNuclearSlope()
{
	return 14.1 * pow(AtomicMass, 2.0 / 3.0);
}

double Material::CalculateRadiationLength()
{
	//Check the required parameters exist!

	//Via the method in the PDG
	double Prefactor = (4 * FineStructureConstant * pow(ElectronRadius, 2) * Avogadro) / GetAtomicMass();
	double Z = GetAtomicNumber();
	double a = FineStructureConstant * Z;
	double a2 = pow(a, 2);

	double L1 = log(184.15 * pow(Z, -1.0 / 3.0));
	double L2 = log(1194 * pow(Z, -2.0 / 3.0));

	double F = a2 * (pow(1 + a2, -1) + 0.20206 - 0.0369 * a2 + 0.0083 * pow(a, 4) - 0.002 * pow(a, 6));

	//remember to put in the factor of the density later
	return 1.0 / (Prefactor * ((Z * Z * (L1 - F)) + (Z * L2)));
}

double Material::CalculateSixtrackTotalNucleusCrossSection()
{
	/*
	 * "Neutron total cross sections on nuclei at Fermilab energies"
	 * P.V.R. Murthy, C.A. Ayre, H.R. Gustafson, L.W. Jones, M.J. Longo
	 * Nuclear Physics B
	 * Volume 92, Issue 3, 16 June 1975, Pages 269–308
	 * http://dx.doi.org/10.1016/0550-3213(75)90182-0
	 * See section 5.3
	 * Note the different power from the inelastic below
	 */

	//No to the above, but the following works roughly:
	return 0.04955 * pow(AtomicMass, 0.77);
}

double Material::CalculateSixtrackInelasticNucleusCrossSection()
{
	/*
	 * "Neutron-nucleus inelastic cross sections from 160 to 375 GeV/c"
	 * T.J. Roberts, H.R. Gustafson, L.W. Jones, M.J. Longo, M.R. Whalley
	 * Nuclear Physics B
	 * Volume 159, Issues 1–2, 5–12 November 1979, Pages 56–66
	 * http://dx.doi.org/10.1016/0550-3213(79)90326-2
	 * Note eq 5 - this should be what the Sixtrack/K2 source (as in the PhD thesis) got the values from.
	 */
	return 0.0412 * pow(AtomicMass, 0.711);
}

double Material::CalculateSixtrackRutherfordCrossSection()
{
	double R = 1.2e-15 * pow(AtomicMass, 1.0 / 3.0);
	double expC = -0.856e-3 * pow(R, 2);
	const double C = pow(PlanckConstantBar * SpeedOfLight / (ElectronCharge * 1e9 * 0.001 * 1e-28), 2);
	const double PiAlpha = 4 * pi * pow(FineStructureConstant, 2);
	return 0;
}

double Material::CalculateSixtrackdEdx(double E)
{
	//Since the numbers in sixtrack make no sense, what can be done here?
	return 1;
}

void Material::Material::SetAtomicNumber(double p)
{
	AtomicNumber = p;
}

void Material::SetAtomicMass(double p)
{
	AtomicMass = p;
}

void Material::SetName(std::string p)
{
	Name = p;
}

void Material::SetSymbol(std::string p)
{
	Symbol = p;
}

void Material::SetConductivity(double p)
{
	Conductivity = p;
}

void Material::SetRadiationLength(double p)
{
	X0 = p;
}
void Material::SetDensity(double p)
{
	Density = p;
}

void Material::SetElectronDensity(double p)
{
	ElectronDensity = p;
}

void Material::SetMeanExcitationEnergy(double p)
{
	MeanExcitationEnergy = p;
}

void Material::SetPlasmaEnergy(double p)
{
	PlasmaEnergy = p;
}

void Material::SetSixtrackTotalNucleusCrossSection(double p)
{
	sigma_pN_total = p;
}

void Material::SetSixtrackInelasticNucleusCrossSection(double p)
{
	sigma_pN_inelastic = p;
}

void Material::SetSixtrackElasticNucleusCrossSection(double p)
{
	sigma_pN_elastic = p;
}

void Material::SetSixtrackRutherfordCrossSection(double p)
{
	sigma_Rutherford = p;
}

void Material::SetSixtrackdEdx(double p)
{
	dEdx = p;
}

void Material::SetSixtrackNuclearSlope(double p)
{
	b_N = p;
}

/*
 * Accessors
 */
double Material::GetAtomicNumber() const
{
	return AtomicNumber;
}

double Material::GetAtomicMass() const
{
	return AtomicMass;
}

std::string Material::GetName() const
{
	return Name;
}

std::string Material::GetSymbol() const
{
	return Symbol;
}

double Material::GetConductivity() const
{
	return Conductivity;
}

double Material::GetRadiationLength() const
{
	return X0;
}

double Material::GetRadiationLengthInM() const
{
	return X0 * 0.001 / Density;
}

double Material::GetDensity() const
{
	return Density;
}

double Material::GetElectronDensity() const
{
	return ElectronDensity;
}

double Material::GetMeanExcitationEnergy() const
{
	return MeanExcitationEnergy;
}

double Material::GetPlasmaEnergy() const
{
	return PlasmaEnergy;
}

double Material::GetSixtrackTotalNucleusCrossSection() const
{
	return sigma_pN_total;
}

double Material::GetSixtrackInelasticNucleusCrossSection() const
{
	return sigma_pN_inelastic;
}

double Material::GetSixtrackElasticNucleusCrossSection() const
{
	return sigma_pN_elastic;
}

double Material::GetSixtrackRutherfordCrossSection() const
{
	return sigma_Rutherford;
}

double Material::GetSixtrackdEdx() const
{
	return dEdx;
}

double Material::GetSixtrackNuclearSlope() const
{
	return b_N;
}

/*
 * Verifies that all the material entries exist and make sense.
 * returns true if the material is good.
 * returns false if the material is bad.
 */
bool Material::VerifyMaterial() const
{
	bool verification = true;
	if(GetName().size() < 1)
	{
		return false;
	}

	if(GetSymbol().size() < 1)
	{
		std::cerr << "Failed to verify Symbol for " << GetName() << ": " << GetSymbol() << std::endl;
		verification = false;
	}
	if(GetAtomicNumber() < 1)
	{
		std::cerr << "Failed to verify AtomicNumber for " << GetName() << ": " << GetAtomicNumber() << std::endl;
		verification = false;
	}
	if(GetAtomicMass() <= 0)
	{
		std::cerr << "Failed to verify AtomicMass for " << GetName() << ": " << GetAtomicMass() << std::endl;
		verification = false;
	}

	if(GetConductivity() <= 0)
	{
		std::cerr << "failed to verify conductivity for " << GetName() << ": " << GetConductivity() << std::endl;
		verification = false;
	}
	if(GetRadiationLength() <= 0)
	{
		std::cerr << "Failed to verify RadiationLength for " << GetName() << ": " << GetRadiationLength() << std::endl;
		verification = false;
	}
	if(GetDensity() <= 0)
	{
		std::cerr << "Failed to verify Density for " << GetName() << ": " << GetDensity() << std::endl;
		verification = false;
	}
	if(GetElectronDensity() <= 0)
	{
		std::cerr << "Failed to verify ElectronDensity for " << GetName() << ": " << GetConductivity() << std::endl;
		verification = false;
	}
	if(GetMeanExcitationEnergy() <= 0)
	{
		std::cerr << "Failed to verify MeanExcitationEnergy for " << GetName() << ": " << GetMeanExcitationEnergy()
				  << std::endl;
		verification = false;
	}
	if(GetPlasmaEnergy() <= 0)
	{
		std::cerr << "Failed to verify PlasmaEnergy for " << GetName() << ": " << GetPlasmaEnergy() << std::endl;
		verification = false;
	}

	/*
	 * Sixtrack parameters
	 */
	if(GetSixtrackTotalNucleusCrossSection() <= 0)
	{
		std::cerr << "Failed to verify SixtrackTotalNucleusCrossSection for " << GetName() << ": "
				  << GetSixtrackTotalNucleusCrossSection() << std::endl;
		verification = false;
	}
	if(GetSixtrackInelasticNucleusCrossSection() <= 0)
	{
		std::cerr << "Failed to verify SixtrackInelasticNucleusCrossSection for " << GetName() << ": "
				  << GetSixtrackInelasticNucleusCrossSection() << std::endl;
		verification = false;
	}
	if(GetSixtrackRutherfordCrossSection() <= 0)
	{
		std::cerr << "Failed to verify SixtrackRutherfordCrossSection for " << GetName() << ": "
				  << GetSixtrackRutherfordCrossSection() << std::endl;
		verification = false;
	}
	if(GetSixtrackNuclearSlope() <= 0)
	{
		std::cerr << "Failed to verify SixtrackNuclearSlope for " << GetName() << ": " << GetSixtrackNuclearSlope()
				  << std::endl;
		verification = false;
	}
	if(GetSixtrackdEdx() <= 0)
	{
		std::cerr << "Failed to verify SixtrackdEdx for " << GetName() << ": " << GetSixtrackdEdx() << std::endl;
		verification = false;
	}

	return verification;
}

bool Material::IsMixture() const
{
	return false;
}

Material* Material::SelectRandomMaterial()
{
	return this;
}
