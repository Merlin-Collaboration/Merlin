#ifndef _COLLIMATOR_DATABASE_HPP_
#define _COLLIMATOR_DATABASE_HPP_
#include <fstream>
#include <string>
#include "Collimators/Material.hpp"
#include "Collimators/Material_Database.hpp"

#include "BeamModel/BeamData.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "RingDynamics/LatticeFunctions.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"

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
	double x_offset;	//Collimator x offset
	double y_offset;	//Collimator y offset
	double j1_tilt;		//Collimator jaw 1 tilt
	double j2_tilt;		//Collimator jaw 2 tilt
	double length;		//Collimator length (m)
	material* Material;	//Collimator material
	double sigma_x;		//Jaw x opening in number of sigmas
	double sigma_y;		//Jaw y opening in number of sigmas
	double beta_x;		//Calculated x beta function at the collimator entrance
	double beta_y;		//Calculated y beta function at the collimator entrance
	double position;	//Length along the lattice, used to calculate the beta functions
};

	collimator* Collimator;
	size_t number_collimators;
	bool use_sigma;

	double Configure_collimators(AcceleratorModel* model, double energy, double emittance_x, double emittance_y, double scale = 1.0e-20);
	void Select_impact_factor(string collimator, double impact);
protected:

	string primary_collimator;		//name of collimator where first impact will occur
	double requested_impact_factor;		//Impact factor in m
	double impact_sigma;			//Impact factor at collimator in number of sigmas

private:

};

#endif
