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

using namespace std;

class MaterialData
{
public:
	map<string, MaterialProperties*> property;
	MaterialData()
	{
	}                     // constructor
	~MaterialData();     // destructor
	void MakeMixture(string name, string s, ...);
	void MakeMixtureByWeight(string name, string s, ...);
	void PrintTable();
};

class StandardMaterialData: public MaterialData
{
public:
	StandardMaterialData();   // constructor fills usual materials
	void UseSixTrackValues();     // there are some differences, goodness knows why
};

#endif /* MATERIALDATA */
