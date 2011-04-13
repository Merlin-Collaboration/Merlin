#include "Collimators/Collimator_Database.hpp"
#include "AcceleratorModel/Apertures/CollimatorAperture.hpp"
#include "Collimators/ResistiveWakePotentials.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

//Read in file into some sensible structure


Collimator_Database::Collimator_Database(string input_file, Material_Database* db, bool sigma) : number_collimators(0),use_sigma(sigma)
{

	//Open file
	ifstream* input = new ifstream(input_file.c_str(), ifstream::in);

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
	requested_impact_factor = 1;
	impact_sigma = 1;
}

double Collimator_Database::Configure_collimators(AcceleratorModel* model, double energy, double emittance_x, double emittance_y, double scale)
{
	//We wish to calculate the lattice functions - alpha/beta
	LatticeFunctionTable* twiss = new LatticeFunctionTable(model,energy);

	//Step size scaling for calculation?
	//Default 1e-16 fails for LHC
	//Use ~1e-20 instead.
	twiss->ScaleBendPathLength(scale);

	//Do the calculation
	twiss->Calculate();

	cout << "Configuring Collimators" << endl;
	vector<Spoiler*> collimators;
	int n_collimators = model->ExtractTypedElements(collimators,"*");
	cout << "Got " << n_collimators << " Collimators" << endl;

	for(vector<Spoiler*>::iterator c = collimators.begin(); c!=collimators.end(); c++)
	{
		//cout <<(*c)->GetName() << "\t" << (*c)->GetLength() << endl;
		for(size_t i=1; i < number_collimators; i++)
		{
			//Time to search for the collimator we are currently using
			if(Collimator[i].name == (*c)->GetName())
			{
				for(int j = 0; j <= twiss->NumberOfRows(); j++)
				{
					if(Collimator[i].position == twiss->Value(0,0,0,j))
					{
						cout << "Hit " << (*c)->GetName() << " at twiss row " << j << " Distance: " << twiss->Value(0,0,0,j);
						double beta_x = twiss->Value(1,1,1,j);		//Beta x
						double beta_y = twiss->Value(3,3,2,j);		//Beta y
						cout << "\t beta x: " << beta_x << "\t beta y: " << beta_y << endl;
						double collimator_aperture_tilt = Collimator[i].tilt;
					if( !(beta_x == 0 || beta_y == 0) )
					{
						double sigma = sqrt( ( beta_x * emittance_x * cos(collimator_aperture_tilt) * cos(collimator_aperture_tilt)) + \
				 				      (beta_y * emittance_y * sin(collimator_aperture_tilt) * sin(collimator_aperture_tilt)) );

						double collimator_aperture_width  = Collimator[i].sigma_x*sigma*2;
						double collimator_aperture_height = Collimator[i].sigma_y*sigma*2;
						if( requested_impact_factor!=1 && (Collimator[i].name == primary_collimator) )
						{
							impact_sigma = ((Collimator[i].sigma_x*sigma+requested_impact_factor)/sigma);
						}

						material* collimator_material = Collimator[i].Material;

						//Create an aperture for the collimator jaws
						CollimatorAperture* app=new CollimatorAperture(collimator_aperture_width,collimator_aperture_height,collimator_aperture_tilt,collimator_material);

						//Set the aperture for collimation
						(*c)->SetAperture(app);


						//Now to set up the resistive wakes
						double conductivity = collimator_material->sigma;
						double aperture_size = collimator_aperture_width;

						//Collimation only will take place on one axis
						if (collimator_aperture_height < collimator_aperture_width)
						{
							aperture_size = collimator_aperture_height;
						} //set to smallest out of height or width

						//Define the resistive wake for the collimator jaws.
						ResistivePotential* resWake =  new ResistivePotential(1,conductivity,0.5*aperture_size,(*c)->GetLength()*meter,"Data/table");

						//Set the Wake potentials for this collimator
						(*c)->SetWakePotentials(resWake);
					}
					else
					{
						cout << "Rejected: " << (*c)->GetName() << endl;
					}

					}
				}
			}
		}
	}
	return impact_sigma;
}

void Collimator_Database::Select_impact_factor(string collimator, double impact)
{
	primary_collimator = collimator;
	requested_impact_factor = impact;		//Impact factor in m
}
