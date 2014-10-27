#ifndef _Material_Database_h_
#define _Material_Database_h_

#include <map>
#include <string>
//#include <vector>
#include "Collimators/Material.h"

using namespace std;

class MaterialDatabase 
{

public:
//Constructor
MaterialDatabase();

//Storage for pointers to material types.
std::map<string,Material*> db;

//Find the material we are interested in
Material* FindMaterial(string symbol);

/*
* Check all the materials in the database are doing something sensible
*/
bool VerifyMaterials();

void DumpMaterialProperties();

private:

};
#endif
