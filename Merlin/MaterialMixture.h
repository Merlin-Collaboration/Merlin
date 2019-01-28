/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _MaterialMixture_h_
#define _MaterialMixture_h_

#include <string>
#include <map>
#include <vector>

#include "Material.h"

/**
 * A material mixture is a mixture of assorted materials, e.g. a metal alloy
 * It contains a map with the component materials.
 * Since it inherits from the base Material class, the same functions can be used in scattering, etc.
 * Some properties will be the weighted mean of the component properties.
 * Some will be discrete values where a material selection must be made at random.
 * To support this the Material functions are virtual, and are overridden here.
 */
class MaterialMixture: public Material
{
public:
	MaterialMixture() :
		Material(), Assembled(false), AssembledByNumber(false), AssembledByMass(false)
	{
	}

	double CalculateElectronDensity();
	double CalculatePlasmaEnergy();
	double CalculateMeanExcitationEnergy();

	double CalculateRadiationLength();
	double CalculateSixtrackdEdx();

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
	std::string GetName() const;
	std::string GetSymbol() const;
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
	 * @return A random element
	 */
	Material* SelectRandomMaterial();

	/**
	 * Returns CurrentElement.
	 * @return CurrentElement
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
	 *
	 * @retval true Corresponds to compounds
	 * @retval false Corresponds to elements
	 */
	virtual bool IsMixture() const;

	/**
	 * Return list of constituent element symbols as strings
	 */
	std::vector<std::pair<std::string, double> > GetConstituentElements();

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
