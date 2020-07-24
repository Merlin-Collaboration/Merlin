/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 *
 *  Created on: 29 Aug 2019
 *      Author: roger
 */

#ifndef MATERIALDATA_H
#define MATERIALDATA_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "MaterialProperties.h"

/**
 * Holds a database of materials and mixtures
 *
 * Each material is stored in a MaterialProperties which can be looked up by name in the map property
 *
 * The StandardMaterialData constructor pre-fills a table with common materials
 */
class MaterialData
{
public:
	std::map<std::string, MaterialProperties*> property;
	MaterialData()
	{
	}                     // constructor
	~MaterialData();     // destructor

	/**
	 * Create a new mixture by combining multiple component with proportions specified by atom count.
	 * @param name Name of mixture
	 * @param s List of components separated by space character
	 * @param proportions List of relative proportions (need not sum to 1)
	 * @param density Density
	 */
	void MakeMixture(std::string name, std::string s, const std::vector<double> &proportions, double density);

	/**
	 * Create a new mixture by combining multiple component with proportions specified by weight.
	 * @param name Name of mixture
	 * @param s List of components separated by space character
	 * @param proportions List of relative proportions by weight (need not sum to 1)
	 * @param density Density
	 */
	void MakeMixtureByWeight(std::string name, std::string s, const std::vector<double> &proportions, double density);

	void MakeMixture(std::string name, std::string s, ...);
	void MakeMixtureByWeight(std::string name, std::string s, ...);
	void PrintTable();
};

/**
 * MaterialData pre-filled with common materials
 */
class StandardMaterialData: public MaterialData
{
public:
	StandardMaterialData();   ///< constructor fills usual materials
	void UseSixTrackValues();     // there are some differences, goodness knows why
};

std::ostream& operator<<(std::ostream& s, MaterialData* M);

#endif /* MATERIALDATA */
