/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _Material_h_
#define _Material_h_

#include <string>

/*
 * Base class for the material definition
 */
class Material
{
public:
	/*
	 * Some of these can be overridden by child classes, e.g. mixtures.
	 * In some cases that doesn't make sense.
	 */

	/**
	 * Overloaded constructor
	 */
	Material()
	{
	}
	Material(std::string name0, std::string sym0, double A0, int AtomicNumber0, double Sigma_E0, double Sigma_I0, double
		Sigma_R0, double dEdx0, double X00, double Density0, double Conductivity0);

	/*
	 * Parameter calculation functions
	 */
	virtual double CalculateElectronDensity();
	virtual double CalculatePlasmaEnergy();
	virtual double CalculateMeanExcitationEnergy();
	virtual double CalculateRadiationLength();

	virtual double CalculateSixtrackNuclearSlope();
	virtual double CalculateSixtrackTotalNucleusCrossSection();
	virtual double CalculateSixtrackInelasticNucleusCrossSection();
	virtual double CalculateSixtrackRutherfordCrossSection();
	virtual double CalculateSixtrackdEdx(double E = 7E12);

	virtual void SetAtomicNumber(double);
	virtual void SetName(std::string);
	virtual void SetSymbol(std::string);
	virtual void SetAtomicMass(double);
	virtual void SetConductivity(double);
	virtual void SetRadiationLength(double);
	virtual void SetDensity(double);
	virtual void SetElectronDensity(double);
	virtual void SetMeanExcitationEnergy(double);
	virtual void SetPlasmaEnergy(double);

	virtual void SetSixtrackTotalNucleusCrossSection(double);
	virtual void SetSixtrackElasticNucleusCrossSection(double);
	virtual void SetSixtrackInelasticNucleusCrossSection(double);
	virtual void SetSixtrackRutherfordCrossSection(double);
	virtual void SetSixtrackdEdx(double);
	virtual void SetSixtrackNuclearSlope(double);

	/*
	 * Define accessors
	 */
	virtual double GetAtomicNumber() const;
	virtual std::string GetName() const;
	virtual std::string GetSymbol() const;
	virtual double GetAtomicMass() const;
	virtual double GetConductivity() const;
	virtual double GetRadiationLength() const;
	virtual double GetRadiationLengthInM() const;
	virtual double GetDensity() const;
	virtual double GetElectronDensity() const;
	virtual double GetMeanExcitationEnergy() const;
	virtual double GetPlasmaEnergy() const;

	virtual double GetSixtrackTotalNucleusCrossSection() const;
	virtual double GetSixtrackInelasticNucleusCrossSection() const;
	virtual double GetSixtrackElasticNucleusCrossSection() const;
	virtual double GetSixtrackRutherfordCrossSection() const;
	virtual double GetSixtrackdEdx() const;
	virtual double GetSixtrackNuclearSlope() const;

	/**
	 * Verifies that all the material entries exist and make sense.
	 * returns true if the material is good.
	 * returns false if the material is bad.
	 */
	virtual bool VerifyMaterial() const;

	/**
	 * Is this a compound material?
	 * true for compounds, false for elements
	 */
	virtual bool IsMixture() const;

	/**
	 * Select a random material for mixtures
	 */
	virtual Material* SelectRandomMaterial();

protected:

	//selected Sixtrack variable names in quotes
	double AtomicNumber;            /// Atomic number (Z)
	double AtomicMass;              /// Atomic mass (A)
	std::string Name;               /// Element/compound name
	std::string Symbol;             /// Elemental symbol: what input file will check for
	double Conductivity;            /// electrical conductivity (sigma) = 1/electrical resisitivity (Ohm*m)e-1 (@ stp)
	double X0;                      /// Radiation Length
	double Density;                 /// Material density, kg/m^3
	double ElectronDensity;         /// Electron density: calculated from other input data
	double MeanExcitationEnergy;
	double PlasmaEnergy;

	//Sixtrack parameters
	double sigma_pN_total;          /// proton nucleus total cross section =sigma_el_pn + sigma_el_pN + sigma_SD_pn + sigma_inel_pN
	double sigma_pN_elastic;        /// proton nucleus elastic cross section
	double sigma_pN_inelastic;      /// proton nucleus inelastic cross section (sigma_inel_pN)
	double sigma_Rutherford;        /// Rutherford scattering cross section (sigma_R)
	double dEdx;                    ///  "dpodx"
	double b_N;                     /// Nuclear slope
};

#endif
