#ifndef _COLLIMATOR_DATABASE_HPP_
#define _COLLIMATOR_DATABASE_HPP_
#include <fstream>
#include <string>
#include "Collimators/Material.hpp"
#include "Collimators/Material_Database.hpp"


using namespace std;
//Collimator database, used to load and store collimator info
class Collimator_Database
{

public:

//Constructor
Collimator_Database(char*, Material_Database*);

struct collimator
{
	string name;		//Collimator name
	float x_gap;		//Collimator x-gap
	float y_gap;		//Collimator y-gap
	float tilt;		//Collimator tilt
	float x_offset;		//Collimator x offset
	float y_offset;		//Collimator y offset
	float j1_tilt;		//Collimator jaw 1 tilt
	float j2_tilt;		//Collimator jaw 2 tilt
	float length;		//Collimator length (m)
	material* Material;	//Collimator material
};

collimator* Collimator;
size_t number_collimators;

protected:

private:

};

#endif
