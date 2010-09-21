#include "Collimators/Collimator_Database.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

//Read in file into some sensible structure


Collimator_Database::Collimator_Database(char* input_file, Material_Database* db, bool sigma) : number_collimators(0),use_sigma(sigma)
{

	//Open file
	ifstream* input = new ifstream(input_file, ifstream::in);

	//Do standard checks
	if(input == NULL || !input->good())
	{
		std::cerr << "Failed to open collimator database: " << input_file << " - Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	//Find number of elements
	string buf;
	if(!use_sigma)
	{
		while(input->good())
		{
			//name	x_gap	y_gap	tilt	x_offset	y_offset	j1_tilt	j2_tilt	length	material
			(*input) >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf;
			number_collimators++;
		}
	}
	else if(use_sigma)
	{
		while(input->good())
		{
			//name	x_sig	y_sig	tilt	material
			(*input) >> buf >> buf >> buf >> buf >> buf;
			number_collimators++;
		}
	}

	std::cout << "found " << number_collimators << " collimators in input file." << std::endl; 

	//Make an appropriately sized array of collimator info
	Collimator = new collimator[number_collimators];

	//Reset the stream
	input->clear();
	input->seekg(0);

	//load into array
	for (unsigned int i=1; i<number_collimators; i++)
	{
		if(!use_sigma)
		{
		//Load up the array of items, put the material name into a string.
		(*input) >> Collimator[i].name >> Collimator[i].x_gap >> Collimator[i].y_gap >> Collimator[i].tilt \
			>> Collimator[i].x_offset >> Collimator[i].y_offset >> Collimator[i].j1_tilt >> Collimator[i].j2_tilt \
			>> Collimator[i].length >> buf;
		}
		else if(use_sigma)
		{
			(*input) >> Collimator[i].name >> Collimator[i].sigma_x >> Collimator[i].sigma_y >> Collimator[i].tilt >> buf;
		}
		//buf contains the material name, must search the material database and check that the appropriate material exists, then adjust the material pointer.
		Collimator[i].Material = db->find_material(buf);
	}


	//Close input file
	input->clear();
	input->close();
	delete input;
	input == NULL;
}
