/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <map>

#include "MaterialMixture.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

#include "RandomNG.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

/*
 * This function adds materials by mass fraction
 * i.e. if 5% of the material by mass is element x, the double here is 0.05, etc.
 */
bool MaterialMixture::AddMaterialByMassFraction(Material* m, double f)
{
	if(AssembledByNumber)
	{
		std::cerr
			<< "Already added a material to this mixture by number fraction. Adding by Mass fraction will break things."
			<< std::endl;
		return false;
	}
	/*
	 * Internally we want to keep track of the number density in order to generate a random nucleus for scattering
	 * The other parameters are constant for the material.
	 */
	std::pair<double, double> p;
	p.first = 0.0;  //Number fraction (currently unknown)
	p.second = f;   //Mass fraction

	std::pair<std::map<Material*, std::pair<double, double> >::iterator, bool> test;
	test = MixtureMap.insert(std::pair<Material*, std::pair<double, double> >(m, p));
	AssembledByMass = true;
	return test.second;
}

/*
 * This function adds materials by number density fraction
 * i.e. if 5% of the component atoms are element m, the double here is 0.05, etc.
 */
bool MaterialMixture::AddMaterialByNumberFraction(Material* m, double f)
{
	if(AssembledByMass)
	{
		std::cerr
			<< "Already added a material to this mixture by mass fraction. Adding by number fraction will break things."
			<< std::endl;
		return false;
	}
	std::pair<double, double> p;
	p.first = f;    //Number fraction
	p.second = 0.0; //Mass fraction (currently unknown)

	std::pair<std::map<Material*, std::pair<double, double> >::iterator, bool> test;
	test = MixtureMap.insert(std::pair<Material*, std::pair<double, double> >(m, p));
	AssembledByNumber = true;
	return test.second;
}

double MaterialMixture::CalculateElectronDensity()
{
	/*
	 * 1: Get weighted mean of the AtomicNumber
	 * 2: Get weighted mean of the AtomicMass
	 */
	double Z = 0.0;
	double A = 0.0;

	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(MaterialIt != MixtureMap.end())
	{
		Z += MaterialIt->second.first * MaterialIt->first->GetAtomicNumber();
		A += MaterialIt->second.first * MaterialIt->first->GetAtomicMass();
		MaterialIt++;
	}
	return Z * Avogadro * GetDensity() * 0.001 / (A * pow(centimeter, 3)); // n_e m^-3 (1e6 conversion from cm^3)
}

//returns the Plasma energy in GeV
double MaterialMixture::CalculatePlasmaEnergy()
{
	return (PlanckConstantBar * sqrt((ElectronDensity * pow(ElectronCharge, 2)) / (ElectronMass
		   * FreeSpacePermittivity))) / ElectronCharge * eV;
}

double MaterialMixture::CalculateMeanExcitationEnergy()
{
	/*
	 * Here we follow the form given in:
	 * "Evaluation of the collision stopping power of elements and compounds for electrons and positrons, Stephen M. Seltzer, Martin J. Berger"
	 * The International Journal of Applied Radiation and Isotopes
	 * Volume 33, Issue 11, November 1982, Pages 1189–1218
	 * http://dx.doi.org/10.1016/0020-708X(82)90244-7
	 * Note eq 10,11
	 * Also note: "The recommendation in Table 6 is to use I-values for these other constituents which are 13% larger than the corresponding I-values for condensed elemental substances in Table 2"
	 * Thus we use I -> 1.13*Ii
	 *
	 * The below is not valid for organic compounds with H,C,N,O,F,Cl.
	 * Look at the values in Table 6 of the above paper instead.
	 */

	double ISum = 0.0;
	double ZA = 0.0;

	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(MaterialIt != MixtureMap.end())
	{
		double A = MaterialIt->first->GetAtomicMass();
		double Z = MaterialIt->first->GetAtomicNumber();
		double I_el = 1.13 * MaterialIt->first->GetMeanExcitationEnergy();
		double w = MaterialIt->second.second;
		ZA += w * Z / A;
		ISum += w * (Z / A) * log(I_el);

		MaterialIt++;
	}
	return exp(ISum / ZA);
}

double MaterialMixture::CalculateSixtrackdEdx()
{
	/*
	 * pdg states for a mixture the following is an approximation:
	 * dE/dx = \sum (w_i * (dE/dx)_i)
	 * where w_i is the mass fraction
	 */

	double dEdx = 0.0;

	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(MaterialIt != MixtureMap.end())
	{
		dEdx += (MaterialIt->second.second * MaterialIt->first->GetSixtrackdEdx());
		MaterialIt++;
	}
	return dEdx;
}

double MaterialMixture::CalculateRadiationLength()
{
	/*
	 * pdg states for a mixture the following is an approximation:
	 * 1/X_0 = \sum (w_i / X_i)
	 * where w_i is the mass fraction
	 */

	double X0 = 0.0;

	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(MaterialIt != MixtureMap.end())
	{
		X0 += (MaterialIt->second.second / MaterialIt->first->GetRadiationLength());
		MaterialIt++;
	}
	return 1.0 / X0;
}

/**
 * Set parameters
 */
void MaterialMixture::SetName(std::string p)
{
	Name = p;
}

void MaterialMixture::SetSymbol(std::string p)
{
	Symbol = p;
}

void MaterialMixture::SetConductivity(double p)
{
	Conductivity = p;
}

void MaterialMixture::SetRadiationLength(double p)
{
	X0 = p;
}

void MaterialMixture::SetDensity(double p)
{
	Density = p;
}

void MaterialMixture::SetElectronDensity(double p)
{
	ElectronDensity = p;
}
/*
   void MaterialMixture::SetElectronCriticalEnergy(double p)
   {
    ElectronCriticalEnergy = p;
   }
 */
void MaterialMixture::SetMeanExcitationEnergy(double p)
{
	MeanExcitationEnergy = p;
}

void MaterialMixture::SetPlasmaEnergy(double p)
{
	PlasmaEnergy = p;
}

void MaterialMixture::SetSixtrackdEdx(double p)
{
	dEdx = p;
}

/**
 * Accessors
 */
std::string MaterialMixture::GetName() const
{
	return Name;
}

std::string MaterialMixture::GetSymbol() const
{
	return Symbol;
}

double MaterialMixture::GetConductivity() const
{
	return Conductivity;
}

double MaterialMixture::GetRadiationLength() const
{
	return X0;
}

double MaterialMixture::GetRadiationLengthInM() const
{
	return X0 * 0.001 / Density;
}

double MaterialMixture::GetDensity() const
{
	return Density;
}

double MaterialMixture::GetElectronDensity() const
{
	return ElectronDensity;
}
/*
   double MaterialMixture::GetElectronCriticalEnergy() const
   {
    return ElectronCriticalEnergy;
   }
 */
double MaterialMixture::GetMeanExcitationEnergy() const
{
	return MeanExcitationEnergy;
}

double MaterialMixture::GetPlasmaEnergy() const
{
	return PlasmaEnergy;
}

double MaterialMixture::GetSixtrackdEdx() const
{
	return dEdx;
}

//Random element
double MaterialMixture::GetAtomicNumber() const
{
	return CurrentMaterial->GetAtomicNumber();
}

//Random element
double MaterialMixture::GetAtomicMass() const
{
	return CurrentMaterial->GetAtomicMass();
}

//Random element
double MaterialMixture::GetSixtrackTotalNucleusCrossSection() const
{
	return CurrentMaterial->GetSixtrackTotalNucleusCrossSection();
}

//Random element
double MaterialMixture::GetSixtrackInelasticNucleusCrossSection() const
{
	return CurrentMaterial->GetSixtrackInelasticNucleusCrossSection();
}

//Random element
double MaterialMixture::GetSixtrackRutherfordCrossSection() const
{
	return CurrentMaterial->GetSixtrackRutherfordCrossSection();
}

//Random element
double MaterialMixture::GetSixtrackNuclearSlope() const
{
	return CurrentMaterial->GetSixtrackNuclearSlope();
}

/*
 * Verifies that all the material entries exist and make sense.
 * returns true if the material is good.
 * returns false if the material is bad.
 */
bool MaterialMixture::VerifyMaterial() const
{
	bool verification = true;

	//Check the mixture related bits
	if(GetName().size() < 1)
	{
		return false;
	}
	if(GetSymbol().size() < 1)
	{
		return false;
	}

	//Then check the component materials
	//Here we need to loop over the component materials
	double MassFraction = 0;
	double NumberFraction = 0;

	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(MaterialIt != MixtureMap.end())
	{
		//Pull the number and mass fractions
		NumberFraction += MaterialIt->second.first;
		MassFraction += MaterialIt->second.second;

		bool test = MaterialIt->first->VerifyMaterial();
		if(test == false)
		{
			verification = false;
			std::cerr << "Material " << GetName() << " subcomponent: "  << MaterialIt->first->GetName()
					  << " failed to verify." << std::endl;
		}
		MaterialIt++;
	}

	if(GetConductivity() <= 0)
	{
		std::cerr << "Failed to verify Conductivity for " << GetName() << ": " << GetConductivity() << std::endl;
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
	//if(GetElectronCriticalEnergy() <=0){std::cerr << "Failed to verify ElectronCriticalEnergy for " << GetName() << ": " << GetElectronCriticalEnergy() << std::endl; verification = false;}
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

	//Verify that the mass and number fraction components sum to 1
	unsigned int oldprec = std::cerr.precision(16);
	double tol = 1e-12;
	if(std::fabs(MassFraction - 1.0) > tol)
	{
		std::cerr << "Invalid mass fraction sum for: " << GetName() << " - " << MassFraction << std::endl;
		verification = false;
	}
	if(std::fabs(NumberFraction - 1.0) > tol)
	{
		std::cerr << "Invalid number fraction sum for: " << GetName() << " - " << NumberFraction << std::endl;
		verification = false;
	}
	std::cerr.precision(oldprec);
	return verification;
}

bool MaterialMixture::Assemble()
{
	std::map<Material*, std::pair<double, double> >::iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	unsigned int count = MixtureMap.size();

	std::vector<double>::iterator FractionVectorIt;
	std::vector<double> FractionVector;
	FractionVector.reserve(count);

	double CurrentFraction, Total = 0.0;
	/*
	 * MassFraction = NumberFraction * AtomicNumber
	 * NumberFraction = MassFraction / AtomicNumber
	 */

	if(AssembledByMass)
	{
		//Set CurrentMaterial with the first element listed
		CurrentMaterial = MaterialIt->first;

		//Want the number fraction here
		//Needs 2 passes
		while(MaterialIt != MixtureMap.end())
		{
			CurrentFraction = MaterialIt->second.second / MaterialIt->first->GetAtomicNumber();
			Total += CurrentFraction;
			FractionVector.push_back(CurrentFraction);
			MaterialIt++;
		}

		MaterialIt = MixtureMap.begin();
		FractionVectorIt = FractionVector.begin();
		while(MaterialIt != MixtureMap.end())
		{
			MaterialIt->second.first = *FractionVectorIt / Total;
			MaterialIt++;
			FractionVectorIt++;
		}
		return true;
	}
	else if(AssembledByNumber)
	{
		//Set CurrentMaterial with the first element listed
		CurrentMaterial = MaterialIt->first;

		while(MaterialIt != MixtureMap.end())
		{
			MaterialIt++;
		}
		return true;
	}
	else
	{
		std::cerr << "Add materials before using Assemble()"  << std::endl;
		return false;
	}
}

Material* MaterialMixture::SelectRandomMaterial()
{
	double x = RandomNG::uniform(0, 1);
	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	MaterialIt = MixtureMap.begin();
	while(x > 0)
	{
		x -= MaterialIt->second.first;
		if(x > 0)
		{
			MaterialIt++;
		}
	}
	//std::cout << "x = " << x << "\t material first =" << MaterialIt->first->GetName() << "\t material second =" << MaterialIt->second.first << std::endl;
	CurrentMaterial = MaterialIt->first;
	return CurrentMaterial;
}

Material* MaterialMixture::GetCurrentMaterial()
{
	return CurrentMaterial;
}

bool MaterialMixture::IsMixture() const
{
	return true;
}

std::vector<std::pair<std::string, double> > MaterialMixture::GetConstituentElements()
{
	std::map<Material*, std::pair<double, double> >::const_iterator MaterialIt;
	std::vector<std::pair<std::string, double> > elements;
	std::pair<std::string, double> test;

	for(MaterialIt = MixtureMap.begin(); MaterialIt != MixtureMap.end(); ++MaterialIt)
	{
		test = make_pair(MaterialIt->first->GetSymbol(), MaterialIt->second.second);
		elements.push_back(test);
	}

	return elements;
}
