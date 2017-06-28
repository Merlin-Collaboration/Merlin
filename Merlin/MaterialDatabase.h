#ifndef _Material_Database_h_
#define _Material_Database_h_

#include <map>
#include <string>
#include <vector>

#include "Material.h"

class MaterialDatabase
{

public:
//Constructor
	MaterialDatabase();

//Storage for pointers to material types.
	std::map<std::string,Material*> db;

//Find the material we are interested in
	Material* FindMaterial(std::string symbol);

	/*
	* Check all the materials in the database are doing something sensible
	*/
	bool VerifyMaterials();

	void DumpMaterialProperties();

private:

};
#endif
