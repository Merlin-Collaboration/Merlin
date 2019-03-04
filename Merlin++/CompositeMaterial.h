/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _CompositeMaterial_h_
#define _CompositeMaterial_h_

#include <string>
#include <map>
#include <vector>

#include "Material.h"

using namespace std;

/**
 * A CompositeMaterial is a composite of assorted materials
 * and contains a map of constituent materials.
 * It is similar to a MaterialMixture however it uses bulk material
 * properties as well as mass fraction weighted values for atomic mass,
 * number, etc, and all cross sections.
 * This makes the CompositeMaterial compatible with CrossSections.
 * Note that scattering will be performed from imaginary composite
 * nuclei rather than the constituent nuclei of the composite, unless
 * this is taken into account by creating CrossSections for each
 * constituent and using the SelectRandomMaterial() function.
 */

class CompositeMaterial: public Material
{
public:
	CompositeMaterial() :
		Material(), Assembled(false), AssembledByNumber(false), AssembledByMass(false)
	{
	}

	double CalculateElectronDensity();
	double CalculatePlasmaEnergy();
	double CalculateMeanExcitationEnergy();

	double CalculateRadiationLength();
	double CalculateSixtrackdEdx();

	/**
	 * This calculates and sets all variables using mass fraction weighting
	 */
	void CalculateAllWeightedVariables();

	// New weighed calculate functions
	double CalculateWeightedA();
	double CalculateWeightedZ();
	double CalculateSixTrackElasticNucleusCrossSection();

	void SetName(std::string);
	void SetSymbol(std::string);
	void SetConductivity(double);
	void SetRadiationLength(double);
	void SetDensity(double);
	void SetElectronDensity(double);
//	void SetElectronCriticalEnergy(double);
	void SetMeanExcitationEnergy(double);
	void SetPlasmaEnergy(double);

	void SetSixtrackdEdx(double);

	//~ void SetAtomicMass(double);
	//~ void SetAtomicNumber(double);
	//~ void SetSixtrackTotalNucleusCrossSection(double);
	//~ void SetSixtrackInelasticNucleusCrossSection(double);
	//~ void SetSixtrackRutherfordCrossSection(double);
	//~ void SetSixtrackNuclearSlope(double);

	// Define accessors
	double GetAtomicNumber() const;
	string GetName() const;
	string GetSymbol() const;
	double GetAtomicMass() const;
	double GetConductivity() const;
	double GetRadiationLength() const;
	double GetRadiationLengthInM() const;
	double GetDensity() const;
	double GetElectronDensity() const;
//	double GetElectronCriticalEnergy() const;
	double GetMeanExcitationEnergy() const;
	double GetPlasmaEnergy() const;

	double GetSixtrackTotalNucleusCrossSection() const;
	double GetSixtrackInelasticNucleusCrossSection() const;
	double GetSixtrackRutherfordCrossSection() const;
	double GetSixtrackdEdx() const;
	double GetSixtrackNuclearSlope() const;

	/**
	 * Check that the material properties make some sort of sense
	 */
	bool VerifyMaterial() const;

	/**
	 * This function adds materials by mass fraction
	 * i.e. if 5% of the material by mass is element m, the double here is 0.05, etc.
	 */
	bool AddMaterialByMassFraction(Material*, double);

	/**
	 * This function adds materials by number density fraction
	 * i.e. if 5% of the component atoms are element m, the double here is 0.05, etc.
	 */
	bool AddMaterialByNumberFraction(Material*, double);

	/**
	 * Returns a random element and sets CurrentElement to this element also.
	 */
	Material* SelectRandomMaterial();

	/**
	 * Returns CurrentElement.
	 */
	Material* GetCurrentMaterial();

	/**
	 * Assembles the material
	 */
	bool Assemble();

	/**
	 * Is this material ready to be used?
	 */
	bool IsAssembled();

	/**
	 * Is this a compound material?
	 * true for compounds, false for elements
	 */
	virtual bool IsMixture() const;

	/**
	 * Return list of constituent element symbols as strings
	 */
	vector<pair<string, double> > GetConstituentElements();

private:

	/**
	 * A map of number density fractions in the material, along with the material pointer.
	 * In the double pair:
	 * first = number fraction
	 * second = mass fraction
	 */
	std::map<Material*, std::pair<double, double> > MixtureMap;

	/**
	 * Current element selected
	 */
	Material* CurrentMaterial;

	bool Assembled;
	bool AssembledByNumber;
	bool AssembledByMass;
};
#endif
