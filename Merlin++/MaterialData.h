/*
 * MaterialData.h
 *
 *  Created on: 29 Aug 2019
 *      Author: roger
 */

#ifndef MATERIALDATA_H
#define MATERIALDATA_H

#include <iostream>
#include <string>
#include <map>
#include "MaterialProperties.h"

class MaterialData
{
public:
	std::map<std::string, MaterialProperties*> property;
	MaterialData()
	{
	}                     // constructor
	~MaterialData();     // destructor
	void MakeMixture(std::string name, std::string s, ...);
	void MakeMixtureByWeight(std::string name, std::string s, ...);
	void PrintTable();
};

class StandardMaterialData: public MaterialData
{
public:
	StandardMaterialData();   // constructor fills usual materials
	void UseSixTrackValues();     // there are some differences, goodness knows why
};

std::ostream& operator<<(std::ostream& s, MaterialData* M);

#endif /* MATERIALDATA */
