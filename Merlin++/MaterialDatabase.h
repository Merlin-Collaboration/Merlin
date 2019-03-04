/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _Material_Database_h_
#define _Material_Database_h_

#include <map>
#include <string>
#include <vector>

#include "Material.h"

class MaterialDatabase
{

public:
	/**
	 * Constructor
	 */
	MaterialDatabase();

	/**
	 * Storage for pointers to material types.
	 */
	std::map<std::string, Material*> db;

	/**
	 * Find the material we are interested in
	 */
	Material* FindMaterial(std::string symbol);

	/**
	 * Check all the materials in the database are doing something sensible
	 */
	bool VerifyMaterials();

	void DumpMaterialProperties();

private:

};
#endif
