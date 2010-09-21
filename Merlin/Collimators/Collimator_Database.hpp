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
Collimator_Database(char*, Material_Database*, bool use_sigma);

struct collimator
{
	string name;		//Collimator name
	double x_gap;		//Collimator x-gap
	double y_gap;		//Collimator y-gap
	double tilt;		//Collimator tilt
	double x_offset;		//Collimator x offset
	double y_offset;		//Collimator y offset
	double j1_tilt;		//Collimator jaw 1 tilt
	double j2_tilt;		//Collimator jaw 2 tilt
	double length;		//Collimator length (m)
	material* Material;	//Collimator material
	double sigma_x;
	double sigma_y;
};

collimator* Collimator;
size_t number_collimators;
bool use_sigma;
protected:

private:

};

#endif
