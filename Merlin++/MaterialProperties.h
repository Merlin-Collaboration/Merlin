/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 *
 *  Created on: 28 Apr 2017
 *      Author: roger
 */

#ifndef MATERIALPROPERTIES_H_
#define MATERIALPROPERTIES_H_

#include <vector>

/**
 * MaterialProperties.h
 *
 *  Contains the physical quantities for a material
 *  Used for collimation and in wakefields
 *
 *  Created on: 16 Aug 2017
 *      Author: roger
 */

#include <iostream>
#include <string>
#include <map>

/**
 * Stores material properties that are used by physics processes.
 *
 * The following are stored as member values and must be specified for all materials
 *
 * Z, A, density, dEdx, sigma_R, sigma_I, sigma_E, sigma_D, sigma_T, X0
 *
 * Additional properties need for specific physics process can be added using
 * SetExtra(), GetExtra() and HaveExtra().
 *
 * For composite materials see Mixture
 */
class MaterialProperties
{
public:
	double Z, A, density, dEdx, sigma_R, sigma_I, sigma_E, sigma_D, sigma_T, X0;
	double lambda;
	std::map<std::string, double> *extra;

	virtual double A_R()
	{
		return A;
	}
	virtual double A_H()
	{
		return A;
	}
	MaterialProperties()
	{
		lambda = A = density = dEdx = sigma_R = sigma_I = sigma_T = Z = X0 = 0;
		extra = nullptr;
	}
	// keep child class happy
	MaterialProperties(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
	void Update();
	MaterialProperties* EnergyScale(double E);
	virtual ~MaterialProperties();
	MaterialProperties(const MaterialProperties&); // copy constructor
	MaterialProperties& operator=(const MaterialProperties&); // copy assignment
	void SetExtra(std::string, ...);
	double GetExtra(std::string);
	bool HaveExtra(std::string);
};

/**
 * Stores a Mixture of 2 or more MaterialProperties.
 */
class Mixture: public MaterialProperties
{
public:
	int N;
	std::vector<double> V_R;
	std::vector<double> V_H;
	std::vector<double> V_A;
	double A_H() override; ///< random nucleus for elastic scattering
	double A_R() override; ///< random nucleus for Rutherford scattering
	Mixture(std::map<std::string, MaterialProperties*> dict, std::vector<std::string> names, std::vector<double>
		proportions, double d);
};

std::ostream& operator<<(std::ostream& o, MaterialProperties& d);

#endif /* MATERIALPROPERTIES_H_ */
