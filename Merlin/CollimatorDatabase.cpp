/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "CollimatorAperture.h"

#include "CollimatorDatabase.h"
#include "ResistiveWakePotentials.h"

#include "PhysicalUnits.h"

#include "RandomNG.h"

using namespace PhysicalUnits;

//Read in file into some sensible structure

CollimatorDatabase::CollimatorDatabase(string input_file, MaterialDatabase* db, bool sigma) :
		number_collimators(0), use_sigma(sigma), logFlag(false), ErrorLogFlag(false), EnableMatchBeamEnvelope(true), EnableMatchReferenceOrbit(
				true), JawFlattnessErrors(false), JawAlignmentErrors(false), EnableResistiveCollimatorWakes(false), AngleError(
				0), PositionError(0)
{

	//Open file
	ifstream* input = new ifstream(input_file.c_str(), ifstream::in);

	//Do standard checks
	if (input == nullptr || !input->good())
	{
		std::cerr << "Failed to open collimator database: " << input_file << " - Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	//Find number of elements
	string buf;
	if (!use_sigma)
	{
		while (input->good())
		{
			//name  x_gap   y_gap   tilt    x_offset    y_offset    j1_tilt j2_tilt length  material
			if ((*input) >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf >> buf)
			{
				number_collimators++;
			}
		}
	}
	else if (use_sigma)
	{
		while (input->good())
		{
			//name  x_sig   y_sig   tilt    material
			if ((*input) >> buf >> buf >> buf >> buf >> buf)
			{
				number_collimators++;
			}
		}
	}

	std::cout << "found " << number_collimators << " collimators in input file." << std::endl;

	//Make an appropriately sized array of collimator info
	CollData = new CollimatorData[number_collimators];

	//Reset the stream
	input->clear();
	input->seekg(0);

	//load into array
	for (unsigned int i = 0; i < number_collimators; i++)
	{
		if (!use_sigma)
		{
			//Load up the array of items, put the material name into a string.
			(*input) >> CollData[i].name >> CollData[i].x_gap >> CollData[i].y_gap >> CollData[i].tilt
					>> CollData[i].x_offset >> CollData[i].y_offset >> CollData[i].j1_tilt >> CollData[i].j2_tilt
					>> CollData[i].length >> buf;
		}
		else if (use_sigma)
		{
			(*input) >> CollData[i].name >> CollData[i].sigma_x >> CollData[i].sigma_y >> CollData[i].tilt >> buf;
		}
		//buf contains the material name, must search the material database and check that the appropriate material exists, then adjust the material pointer.
		CollData[i].JawMaterial = db->FindMaterial(buf);
	}

	//Close input file
	input->clear();
	input->close();
	delete input;
	RequestedImpactFactor = 1;
	ImpactSigma = 1;
}

void CollimatorDatabase::ConfigureCollimators(AcceleratorModel* model)
{
	vector<Collimator*> Collimators;
	size_t n_collimators = model->ExtractTypedElements(Collimators, "*");
	if (logFlag)
	{
		log->setf(ios::left);
		//*log << "Configuring Collimators" << endl;
		*log << "#Got " << n_collimators << " Collimators" << endl;
		*log << "#" << std::setw(20) << "name" << std::setw(7) << "row" << std::setw(14) << std::setw(15) << "Distance"
				<< std::setw(15) << "beta x "\
 << std::setw(15) << "beta y" << std::setw(15) << "cent x" << std::setw(15)
				<< "cent y" << std::setw(15) << "x half gap" << std::setw(15) << "y half gap" << endl;
	}
	map<string, Collimator*> CollimatorMap;
	map<string, Collimator*>::iterator CMapit;

	for (vector<Collimator*>::iterator c = Collimators.begin(); c != Collimators.end(); c++)
	{
		CollimatorMap.insert(pair<string, Collimator*>((*c)->GetName(), (*c)));
	}

	for (size_t n = 0; n < number_collimators; n++)
	{
		CMapit = CollimatorMap.find(CollData[n].name);
		if (CMapit != CollimatorMap.end())
		{
			//cout << "Found: " << (CMapit->second)->GetQualifiedName() << endl;
			Material* collimator_material = CollData[n].JawMaterial;
			//Create an aperture for the collimator jaws
			CollimatorAperture* app = new CollimatorAperture(CollData[n].x_gap, CollData[n].y_gap, CollData[n].tilt,
					(CMapit->second)->GetLength(), 0, 0);

			//if(EnableMatchBeamEnvelope)
			//if(true)
			{
				app->SetExitWidth(CollData[n].x_gap);      //Horizontal
				app->SetExitHeight(CollData[n].y_gap);    //Vertical
				app->SetExitXOffset(0); //Horizontal
				app->SetExitYOffset(0); //Vertical
			}
			(CMapit->second)->SetAperture(app);
			(CMapit->second)->SetMaterial(collimator_material);
		}
	}
}
double CollimatorDatabase::ConfigureCollimators(AcceleratorModel* model, double emittance_x, double emittance_y,
		LatticeFunctionTable* twiss)
{
	//std::cout << "point1" << std::endl;
	vector<Collimator*> Collimators;
	size_t n_collimators = model->ExtractTypedElements(Collimators, "*");
	//std::cout << "point2" << std::endl;
	if (logFlag)
	{
		log->setf(ios::left);
		//*log << "Configuring Collimators" << endl;
		*log << "#Got " << n_collimators << " Collimators" << endl;
		*log << "#" << std::setw(20) << "name" << std::setw(7) << "row" << std::setw(14) << std::setw(15) << "Distance"
				<< std::setw(15) << "beta x "\
 << std::setw(15) << "beta y" << std::setw(15) << "cent x" << std::setw(15)
				<< "cent y" << std::setw(15) << "x half gap" << std::setw(15) << "y half gap" << endl;
	}
	map<string, Collimator*> CollimatorMap;
	map<string, Collimator*>::iterator CMapit;
	//std::cout << "point3" << std::endl;
	for (vector<Collimator*>::iterator c = Collimators.begin(); c != Collimators.end(); c++)
	{
		CollimatorMap.insert(pair<string, Collimator*>((*c)->GetName(), (*c)));
		//cout << (*c)->GetName() << endl;
	}
	//std::cout << "point4" << std::endl;
	//for(vector<Collimator*>::iterator c = collimators.begin(); c!=collimators.end(); c++)
	//{
	//cout <<(*c)->GetName() << "\t" << (*c)->GetLength() << endl;
	for (size_t i = 0; i < number_collimators; i++)
	{
		//Time to search for the collimator we are currently using
		//if(CollData[i].name == (*c)->GetName())
		//{
		CMapit = CollimatorMap.find(CollData[i].name);
		if (CMapit != CollimatorMap.end())
		{
			for (int j = 0; j <= twiss->NumberOfRows(); j++)
			{
				//                  if(CollData[i].position == twiss->Value(0,0,0,j))
				//std::cout << "point5 - " << twiss->NumberOfRows() << std::endl;
				//std::cout << (CMapit->second)->GetComponentLatticePosition() << "\t-\t" << twiss->Value(0,0,0,j) << std::endl;
				if ((CMapit->second)->GetComponentLatticePosition() == twiss->Value(0, 0, 0, j))
				{
					if (logFlag)
					{
						*log << std::setw(20) << (CMapit->second)->GetName() << std::setw(7) << j << std::setw(14)
								<< twiss->Value(0, 0, 0, j);
					}
					//*log << "Found " << std::setw(25) << (*c)->GetName() << " at twiss row " << std::setw(6) << j << " Distance: " << std::setw(7) << twiss->Value(0,0,0,j);
					//std::setw(20) << (CMapit->second)->GetName() << std::setw(7) << j << std::setw(14) << twiss->Value(0,0,0,j);
					//*log << "Found " << std::setw(25) << (*c)->GetName() << " at twiss row " << std::setw(6) << j << " Distance: " << std::setw(7) << twiss->Value(0,0,0,j);
					if (ErrorLogFlag)
					{
						*ErrorLog << std::setw(20) << (CMapit->second)->GetName() << std::setw(7) << j;
					}

					double beta_x = twiss->Value(1, 1, 1, j);      //Beta x
					double beta_y = twiss->Value(3, 3, 2, j);      //Beta y
					//std::cout << "beta_x" << beta_x << std::endl;
					//Added x and y orbit parameters to center collimators on where the beam actually goes.
					double x_orbit = twiss->Value(1, 0, 0, j);
					double y_orbit = twiss->Value(3, 0, 0, j);
					double collimator_aperture_tilt = CollData[i].tilt;

					//~ cout << "CollDB: TILT = " << CollData[i].tilt << endl;
					//~ cout << "CollDB: TILT = " << collimator_aperture_tilt << endl;

					//New - we also want to align the collimator to the reference orbit and the beta function
					double beta_x_exit = twiss->Value(1, 1, 1, j + 1);       //Beta x
					double beta_y_exit = twiss->Value(3, 3, 2, j + 1);       //Beta y

					//Added x and y orbit parameters to center collimators on where the beam actually goes.
					double x_orbit_exit = twiss->Value(1, 0, 0, j + 1);
					double y_orbit_exit = twiss->Value(3, 0, 0, j + 1);
					//std::cout << "x_orbit_exit" << x_orbit_exit << std::endl;
					if (logFlag)
					{
						*log << std::setw(15) << beta_x << std::setw(15) << beta_y << std::setw(15) << x_orbit
								<< std::setw(15) << y_orbit;
					}
					//std:: cout << std::setw(15) << beta_x << std::setw(15)<< beta_y << std::setw(15)<< x_orbit << std::setw(15) << y_orbit;
					if (!(beta_x == 0 || beta_y == 0))
					{
						double sigma_entrance = sqrt(
								(beta_x * emittance_x * cos(collimator_aperture_tilt) * cos(collimator_aperture_tilt))
										+ (beta_y * emittance_y * sin(collimator_aperture_tilt)
												* sin(collimator_aperture_tilt)));

						double collimator_aperture_width_entrance = CollData[i].sigma_x * sigma_entrance * 2;
						double collimator_aperture_height_entrance = CollData[i].sigma_y * sigma_entrance * 2;

						//And the exit parameters
						double sigma_exit = sqrt(
								(beta_x_exit * emittance_x * cos(collimator_aperture_tilt)
										* cos(collimator_aperture_tilt))
										+ (beta_y_exit * emittance_y * sin(collimator_aperture_tilt)
												* sin(collimator_aperture_tilt)));

						double collimator_aperture_width_exit = CollData[i].sigma_x * sigma_exit * 2;
						double collimator_aperture_height_exit = CollData[i].sigma_y * sigma_exit * 2;

						//Get the length
						double length = (CMapit->second)->GetLength();
						//std::cout << "length - " << length << std::endl;
						if (RequestedImpactFactor != 1 && (CollData[i].name == PrimaryCollimator))
						{
							ImpactSigma = ((CollData[i].sigma_x * sigma_entrance + RequestedImpactFactor)
									/ sigma_entrance);
							cout << "CollimatorDatabase::ConfigureCollimators : Beta_x = " << beta_x << " Sigma_x = "
									<< sqrt(beta_x * emittance_x) << endl;
							cout << "CollimatorDatabase::ConfigureCollimators : Beta_y = " << beta_y << " Sigma_y = "
									<< sqrt(beta_y * emittance_y) << endl;
							cout << "CollimatorDatabase::ConfigureCollimators : orbit_x = " << x_orbit << " orbit_y = "
									<< y_orbit << endl;

						}

						Material* collimator_material = CollData[i].JawMaterial;

						(CMapit->second)->SetCollID(i + 1);

						FlukaData* fluka_data = new FlukaData;
						fluka_data->id_coll = (CMapit->second)->GetCollID();
						fluka_data->name = CollData[i].name;
						fluka_data->position = (CMapit->second)->GetComponentLatticePosition();
						fluka_data->angle = CollData[i].tilt;
						fluka_data->beta_x = beta_x;
						fluka_data->beta_y = beta_y;
						fluka_data->half_gap = CollData[i].sigma_x * sigma_entrance;
						fluka_data->material = collimator_material->GetSymbol();
						fluka_data->length = length;
						fluka_data->sig_x = sqrt(emittance_x * beta_x);
						fluka_data->sig_y = sqrt(emittance_y * beta_y);
						fluka_data->j1_tilt = 0.;
						fluka_data->j2_tilt = 0.;
						fluka_data->n_sig = CollData[i].sigma_x;

						StoredFlukaData.push_back(fluka_data);

						//std::cout << "EnableMatchBeamEnvelope - " << EnableMatchBeamEnvelope << "\t" << JawFlattnessErrors << std::endl;
						//Create an aperture for the collimator jaws
						if (EnableMatchBeamEnvelope && !JawFlattnessErrors && !JawAlignmentErrors)
						{
							if (logFlag)
							{
								*log << std::setw(15) << collimator_aperture_width_entrance / 2.0 << std::setw(15)
										<< collimator_aperture_height_entrance / 2.0 << endl;
							}

							CollimatorAperture* app = new CollimatorAperture(collimator_aperture_width_entrance,
									collimator_aperture_height_entrance, collimator_aperture_tilt, length, x_orbit,
									y_orbit);

							app->SetExitWidth(collimator_aperture_width_exit);  //Horizontal
							app->SetExitHeight(collimator_aperture_height_exit);    //Vertical
							app->SetExitXOffset(x_orbit_exit);  //Horizontal
							app->SetExitYOffset(y_orbit_exit);  //Vertical
							//Set the aperture for collimation
							(CMapit->second)->SetAperture(app);
							(CMapit->second)->SetMaterial(collimator_material);
						}
						else if (!EnableMatchBeamEnvelope && !JawFlattnessErrors && !JawAlignmentErrors)
						{

							//We will want to calculate the position of the left and right jaws
							//First
							double gap_x_entrance = collimator_aperture_width_entrance / 2;
							double x_rot = (x_orbit * cos(-collimator_aperture_tilt))
									- (y_orbit * sin(-collimator_aperture_tilt));
							x_rot = x_orbit;
							double x1 = x_rot + gap_x_entrance;
							double x2 = x_rot - gap_x_entrance;

							//Second
							double gap_y_entrance = collimator_aperture_height_entrance / 2;
							double y_rot = (x_orbit * sin(-collimator_aperture_tilt))
									+ (y_orbit * cos(-collimator_aperture_tilt));
							y_rot = y_orbit;
							double y1 = y_rot + gap_y_entrance;
							double y2 = y_rot - gap_y_entrance;

							//EXIT
							//First
							double gap_x_exit = collimator_aperture_width_exit / 2;
							double x_rot_exit = (x_orbit_exit * cos(-collimator_aperture_tilt))
									- (y_orbit_exit * sin(-collimator_aperture_tilt));
							x_rot_exit = x_orbit_exit;
							double xx1 = x_rot_exit + gap_x_exit;
							double xx2 = x_rot_exit - gap_x_exit;

							//Second
							double gap_y_exit = collimator_aperture_height_exit / 2;
							double y_rot_exit = (x_orbit_exit * sin(-collimator_aperture_tilt))
									+ (y_orbit_exit * cos(-collimator_aperture_tilt));
							y_rot_exit = y_orbit_exit;
							double yy1 = y_rot_exit + gap_y_exit;
							double yy2 = y_rot_exit - gap_y_exit;

							//cout << (CMapit->second)->GetName() << "\t" << x1 << "\t" << x2 << "\t" << y1 << "\t" << y2 << "\t" <<gap_x_entrance << "\t" <<gap_y_entrance << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<gap_x_entrance << "\t" <<gap_y_entrance << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<(x1-xx1)*2*1000*1000 << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<x2*2 << "\t" << xx2*2<< endl;

							double xj1 = max(x1, xx1);   //+ve values
							double xj2 = min(x2, xx2);   //-ve values

							double yj1 = max(y1, yy1);   //+ve values
							double yj2 = min(y2, yy2);   //-ve values

							double x_size = (xj1 - xj2);
							double y_size = (yj1 - yj2);

							double x_pos = (xj1 + xj2) / 2;
							double y_pos = (yj1 + yj2) / 2;

							//cout << (CMapit->second)->GetName() << "\t" << x_size <<"\t" <<  y_size << "\t" <<x_pos << "\t" << y_pos<< endl;
							//Left jaw, right jaw?

							//*log << std::setw(15) << xj1 << std::setw(15)<< xj2 <<std::setw(15)<< x_pos << endl;
							if (CollData[i].name == "TCDQA.A4R6.B1")
							{
								OneSidedUnalignedCollimatorAperture* app = new OneSidedUnalignedCollimatorAperture(
										x_size, y_size, collimator_aperture_tilt, length, x_pos, y_pos);

								//Set the aperture for collimation
								(CMapit->second)->SetAperture(app);
								(CMapit->second)->SetMaterial(collimator_material);
								//cout << "collimator: " << CollData[i].name << "has aperture" << Collimators[i]->GetName() << endl;
							}
							else if (CollData[i].name == "TCDQA.B4R6.B1")
							{
								OneSidedUnalignedCollimatorAperture* app = new OneSidedUnalignedCollimatorAperture(
										x_size, y_size, collimator_aperture_tilt, length, x_pos, y_pos);

								//Set the aperture for collimation
								(CMapit->second)->SetAperture(app);
								(CMapit->second)->SetMaterial(collimator_material);
								//cout << "collimator: " << CollData[i].name << "has aperture" << Collimators[i]->GetName() << endl;
							}
							else if (CollData[i].name == "TCDQA.C4R6.B1")
							{
								OneSidedUnalignedCollimatorAperture* app = new OneSidedUnalignedCollimatorAperture(
										x_size, y_size, collimator_aperture_tilt, length, x_pos, y_pos);

								//Set the aperture for collimation
								(CMapit->second)->SetAperture(app);
								(CMapit->second)->SetMaterial(collimator_material);
								//cout << "collimator: " << CollData[i].name << "has aperture" << Collimators[i]->GetName() << endl;
							}
							else
							{
								UnalignedCollimatorAperture* app = new UnalignedCollimatorAperture(x_size, y_size,
										collimator_aperture_tilt, length, x_pos, y_pos);

								//Set the aperture for collimation
								(CMapit->second)->SetAperture(app);
								(CMapit->second)->SetMaterial(collimator_material);
								//cout << "COLLIMATOR: " << CollData[i].name << "has aperture" << Collimators[i]->GetName() << endl;
							}
							if (logFlag)
							{
								*log << std::setw(15) << collimator_aperture_width_entrance / 2.0 << std::setw(15)
										<< collimator_aperture_height_entrance / 2.0 << endl;
							}
							//std::cout << std::setw(15) << collimator_aperture_width_entrance/2.0 << std::setw(15)<< collimator_aperture_height_entrance/2.0 << endl;
						}
						else if (!EnableMatchBeamEnvelope && !JawFlattnessErrors && JawAlignmentErrors)
						{

							//std::cout << "point6" << std::endl;
#ifdef ENABLE_MPI
							//Lets start with Jaw 1
							//Random x,y
							int MPI_RANK = MPI::COMM_WORLD.Get_rank();
							int MPI_SIZE = MPI::COMM_WORLD.Get_size();

							double xOffsetError1;
							double yOffsetError1;
							double xAngleError1;
							double yAngleError1;
							double xOffsetError2;
							double yOffsetError2;
							double xAngleError2;
							double yAngleError2;
							double wholeOffsetError;

							if(MPI_RANK == 0)
							{

								xOffsetError1 = RandomNG::normal(0,PositionError,3);
								yOffsetError1 = 0;					//RandomNG::uniform(-PositionError,PositionError);

								//Random theta1, theta2 - small angle approx
								xAngleError1 = length * RandomNG::normal(0,AngleError,3);
								yAngleError1 = 0;//length * RandomNG::uniform(-AngleError,AngleError);

								//Jaw 2
								//Random x,y
								xOffsetError2 = RandomNG::normal(0,PositionError,3);
								yOffsetError2 = 0;//RandomNG::uniform(-PositionError,PositionError);

								//Random theta1, theta2 - small angle approx
								xAngleError2 = length *RandomNG::normal(0,AngleError,3);
								yAngleError2 = 0;//length * RandomNG::uniform(-AngleError,AngleError);
								//cout << "xOffset1\t" << xOffsetError1 << "\t" << "xOffset2\t" << xOffsetError2 <<"\t xAngleError1\t" << xAngleError1 << "\t" << "xAngleError2\t" << xAngleError2 << endl;
								//pointless sync point
								MPI::COMM_WORLD.Barrier();
								double ErrorArray[8];
								ErrorArray[0] = xOffsetError1;
								ErrorArray[1] = yOffsetError1;
								ErrorArray[2] = xAngleError1;
								ErrorArray[3] = yAngleError1;
								ErrorArray[4] = xOffsetError2;
								ErrorArray[5] = yOffsetError2;
								ErrorArray[6] = xAngleError2;
								ErrorArray[7] = yAngleError2;

								//Send to nodes
								for(int n = 1; n < MPI_SIZE; n++)
								{
									MPI::COMM_WORLD.Send(&ErrorArray, 8, MPI::DOUBLE, n, 1);
								}
							}
							else
							{
								MPI::COMM_WORLD.Barrier();
								//Make a buffer
								double ErrorArray[8];

								//Recv new momentum
								MPI::COMM_WORLD.Recv(&ErrorArray, 8, MPI::DOUBLE, 0, 1);
								xOffsetError1 = ErrorArray[0];
								yOffsetError1 = ErrorArray[1];
								xAngleError1 = ErrorArray[2];
								yAngleError1 = ErrorArray[3];
								xOffsetError2 = ErrorArray[4];
								yOffsetError2 = ErrorArray[5];
								xAngleError2 = ErrorArray[6];
								yAngleError2 = ErrorArray[7];
							}
#endif

#ifndef ENABLE_MPI
							double wholeOffsetError;
							double xOffsetError1 = RandomNG::normal(0, PositionError, 3);
							double yOffsetError1 = 0;				//RandomNG::uniform(-PositionError,PositionError);

							//Random theta1, theta2 - small angle approx
							double xAngleError1 = length * RandomNG::normal(0, AngleError, 3);
							double yAngleError1 = 0;			//length * RandomNG::uniform(-AngleError,AngleError);

							//Jaw 2
							//Random x,y
							double xOffsetError2 = RandomNG::normal(0, PositionError, 3);
							double yOffsetError2 = 0;				//RandomNG::uniform(-PositionError,PositionError);
							//cout << "xOffsetError1" << "\t" << xOffsetError1 << "\t" << "xOffsetError2" << "\t" << xOffsetError2 << endl;

							//Random theta1, theta2 - small angle approx
							double xAngleError2 = length * RandomNG::normal(0, AngleError, 3);
							double yAngleError2 = 0;			//length * RandomNG::uniform(-AngleError,AngleError);
#endif

							//ENTRANCE
							//We will want to calculate the position of the left and right jaws
							//First
							double gap_x_entrance = collimator_aperture_width_entrance / 2;
							double x_rot = (x_orbit * cos(-collimator_aperture_tilt))
									- (y_orbit * sin(-collimator_aperture_tilt));
							x_rot = x_orbit;
							wholeOffsetError = RandomNG::normal(0, 2.5 * nanometer, 3);
							cout << "wholeOffsetError" << "\t" << wholeOffsetError << endl;
							double x1 = x_rot + gap_x_entrance + xOffsetError1 + wholeOffsetError;
							double x2 = x_rot - gap_x_entrance + xOffsetError2 + wholeOffsetError;

							//Second
							double gap_y_entrance = collimator_aperture_height_entrance / 2;
							double y_rot = (x_orbit * sin(-collimator_aperture_tilt))
									+ (y_orbit * cos(-collimator_aperture_tilt));
							y_rot = y_orbit;

							double y1 = y_rot + gap_y_entrance + yOffsetError1;
							double y2 = y_rot - gap_y_entrance + yOffsetError2;

							//EXIT
							//First
							//double gap_x_exit = collimator_aperture_width_exit/2;
							double gap_x_exit = collimator_aperture_width_entrance / 2;
							double x_rot_exit = (x_orbit_exit * cos(-collimator_aperture_tilt))
									- (y_orbit_exit * sin(-collimator_aperture_tilt));
							x_rot_exit = x_orbit_exit;

							//double xx1 = x_rot_exit + gap_x_exit + xOffsetError1 + xAngleError1;
							//double xx2 = x_rot_exit - gap_x_exit + xOffsetError2 + xAngleError2;
							//cout << "xx1" << "\t" << xx1 << "\t" << "xx2" << "\t" << xx2  << endl;

							double xx1 = x_rot_exit + gap_x_exit + xOffsetError1 + xAngleError1 + wholeOffsetError;
							double xx2 = x_rot_exit - gap_x_exit + xOffsetError2 + xAngleError2 + wholeOffsetError;

							//Second
							double gap_y_exit = collimator_aperture_height_exit / 2;
							double y_rot_exit = (x_orbit_exit * sin(-collimator_aperture_tilt))
									+ (y_orbit_exit * cos(-collimator_aperture_tilt));
							y_rot_exit = y_orbit_exit;

							//double yy1 = yy_rot + gap_y_exit;
							//double yy2 = yy_rot - gap_y_exit;

							//double yy1 = y_rot_exit + gap_y_exit + yOffsetError1 + yAngleError1;
							//double yy2 = y_rot_exit - gap_y_exit + yOffsetError2 + yAngleError2;
							double yy1 = y_rot_exit + gap_y_exit + yOffsetError1;
							double yy2 = y_rot_exit - gap_y_exit + yOffsetError2;

							//cout << (CMapit->second)->GetName() << "\t" << x_rot << "\t" << x_rot_exit << "\t" << gap_x_entrance << "\t" << gap_x_exit << "\t" << collimator_aperture_width_entrance << "\t" <<collimator_aperture_width_exit << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1 << "\t" << x2 << "\t" << xx1 << "\t" << xx2 << "\t" <<gap_x_entrance << "\t" <<gap_y_entrance << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<gap_x_entrance << "\t" <<gap_y_entrance << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<(x1-xx1)*2*1000*1000 << endl;
							//cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<x2*2 << "\t" << xx2*2<< endl;

							//                          double xj1 = max(x1,xx1);   //+ve values
							//                          double xj2 = min(x2,xx2);   //-ve values

							//                          double yj1 = max(y1,yy1);   //+ve values
							//                          double yj2 = min(y2,yy2);   //-ve values

							double x_size_entrance = (x1 - x2);
							double y_size_entrance = (y1 - y2);

							double x_pos_entrance = (x1 + x2) / 2;
							double y_pos_entrance = (y1 + y2) / 2;

							double x_size_exit = (xx1 - xx2);
							double y_size_exit = (yy1 - yy2);

							double x_pos_exit = (xx1 + xx2) / 2;
							double y_pos_exit = (yy1 + yy2) / 2;

							//cout << "point star-" <<  (CMapit->second)->GetName() << endl;
							//Left jaw, right jaw?
							if (logFlag)
							{
								*log << std::setw(15) << collimator_aperture_width_entrance / 2.0 << std::setw(15)
										<< collimator_aperture_height_entrance / 2.0 << endl;
							}
							//*log << std::setw(15) << xj1 << std::setw(15)<< xj2 <<std::setw(15)<< x_pos << endl;

							if (ErrorLogFlag)
								*ErrorLog << std::setw(15) << xOffsetError1 / micrometer << std::setw(15)
										<< yOffsetError1 / micrometer << std::setw(15) << xAngleError1 / microradian
										<< std::setw(15) << yAngleError1 / microradian << std::setw(15)
										<< xOffsetError2 / micrometer << std::setw(15) << yOffsetError2 / micrometer
										<< std::setw(15) << xAngleError2 / microradian << std::setw(15)
										<< yAngleError2 / microradian << endl;

							//UnalignedCollimatorAperture* app=new UnalignedCollimatorAperture(collimator_aperture_width_entrance,collimator_aperture_height_entrance,

							//x_size += xOffsetError1;
							//y_size += yOffsetError1;

							CollimatorAperture* app = new CollimatorAperture(x_size_entrance, y_size_entrance,
									collimator_aperture_tilt, length, x_pos_entrance, y_pos_entrance);

							app->SetExitXOffset(x_pos_exit);    //Horizontal
							app->SetExitYOffset(y_pos_exit);    //Vertical

							app->SetExitWidth(x_size_exit); //Horizontal
							app->SetExitHeight(y_size_exit);    //Vertical

							//Set the aperture for collimation
							(CMapit->second)->SetAperture(app);
							(CMapit->second)->SetMaterial(collimator_material);
						}
						else
						{
							std::cerr << "Unsupported collimator configuration in CollimatorDatabase.cpp" << std::endl;
							exit(1);
						}

						//std::cout << "point7" << std::endl;
						//Now to set up the resistive wakes
						double conductivity = collimator_material->GetConductivity();
						double aperture_size = collimator_aperture_width_entrance;

						//Collimation only will take place on one axis
						if (collimator_aperture_height_entrance < collimator_aperture_width_entrance)
						{
							aperture_size = collimator_aperture_height_entrance;
						} //set to smallest out of height or width

						if (EnableResistiveCollimatorWakes)
						{
							//Define the resistive wake for the collimator jaws.
							ResistivePotential* resWake = new ResistivePotential(1, conductivity, 0.5 * aperture_size,
									length * meter, "Data/table");

							//Set the Wake potentials for this collimator
							(CMapit->second)->SetWakePotentials(resWake);
						}
					}
					else
					{
						if (logFlag)
						{
							*log << "Rejected: " << (CMapit->second)->GetName() << endl;
						}
						//std::cout << "Rejected "  << std::endl;
					}

				}
			}
		}
		//}
	}
	if (logFlag)
	{
		*log << "Found " << PrimaryCollimator << "\tImpact sigma: " << ImpactSigma << endl;
	}
	return ImpactSigma;
}
/*
 void CollimatorDatabase::SetOneSideTCDQA(AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
 {

 vector<Collimator*> Collimators;
 size_t n_collimators = model->ExtractTypedElements(Collimators,"TCDQA*");
 map<string,Collimator*> CollimatorMap;
 map<string,Collimator*>::iterator CMapit;
 for(vector<Collimator*>::iterator ch = Collimators.begin(); ch!=Collimators.end(); ch++)
 {
 CollimatorMap.insert(pair<string,Collimator*>((*ch)->GetName(),(*ch)));
 }
 //cout <<(*ch)->GetName() << "\t" << (*ch)->GetLength() << endl;
 for(size_t i=0; i < number_collimators; i++)
 {
 //Time to search for the collimator we are currently using
 //      if(CollData[i].name == (*c)->GetName())
 //      {
 CMapit = CollimatorMap.find(CollData[i].name);
 if(CMapit != CollimatorMap.end())
 {
 for(int j = 0; j <= twiss->NumberOfRows(); j++)
 {
 //                  if(CollData[i].position == twiss->Value(0,0,0,j))
 if((CMapit->second)->GetComponentLatticePosition() == twiss->Value(0,0,0,j))
 {
 double beta_x = twiss->Value(1,1,1,j);      //Beta x
 double beta_y = twiss->Value(3,3,2,j);      //Beta y

 //Added x and y orbit parameters to center collimators on where the beam actually goes.
 double x_orbit = twiss->Value(1,0,0,j);
 double y_orbit = twiss->Value(3,0,0,j);
 double collimator_aperture_tilt = Collimator[i].tilt;

 //New - we also want to align the collimator to the reference orbit and the beta function
 double beta_x_exit = twiss->Value(1,1,1,j+1);       //Beta x
 double beta_y_exit = twiss->Value(3,3,2,j+1);       //Beta y

 //Added x and y orbit parameters to center collimators on where the beam actually goes.
 double x_orbit_exit = twiss->Value(1,0,0,j+1);
 double y_orbit_exit = twiss->Value(3,0,0,j+1);

 if(logFlag)
 *log << std::setw(15) << beta_x << std::setw(15)<< beta_y << std::setw(15)<< x_orbit << std::setw(15) << y_orbit;

 if( !(beta_x == 0 || beta_y == 0) )
 {
 double sigma_entrance = sqrt( ( beta_x * emittance_x * cos(collimator_aperture_tilt) * cos(collimator_aperture_tilt)) + \
                                      (beta_y * emittance_y * sin(collimator_aperture_tilt) * sin(collimator_aperture_tilt)) );

 double collimator_aperture_width_entrance  = CollData[i].sigma_x*sigma_entrance*2;
 double collimator_aperture_height_entrance = CollData[i].sigma_y*sigma_entrance*2;

 //And the exit parameters
 double sigma_exit = sqrt( ( beta_x_exit * emittance_x * cos(collimator_aperture_tilt) * cos(collimator_aperture_tilt)) + \
                                      (beta_y_exit * emittance_y * sin(collimator_aperture_tilt) * sin(collimator_aperture_tilt)) );

 double collimator_aperture_width_exit  = CollData[i].sigma_x*sigma_exit*2;
 double collimator_aperture_height_exit = CollData[i].sigma_y*sigma_exit*2;

 //Get the length
 double length = (CMapit->second)->GetLength();

 Material* collimator_material = CollData[i].JawMaterial;


 #ifdef ENABLE_MPI
 //Lets start with Jaw 1
 //Random x,y
 int MPI_RANK = MPI::COMM_WORLD.Get_rank();
 int MPI_SIZE = MPI::COMM_WORLD.Get_size();

 double xOffsetError1;
 double yOffsetError1;
 double xAngleError1;
 double yAngleError1;
 double xOffsetError2;
 double yOffsetError2;
 double xAngleError2;
 double yAngleError2;
 double wholeOffsetDisp(0.1);

 #endif

 #ifndef ENABLE_MPI
 double wholeOffsetError;
 double xOffsetError1 = 0;
 double yOffsetError1 = 0;//RandomNG::uniform(-PositionError,PositionError);

 //Random theta1, theta2 - small angle approx
 double xAngleError1 = 0;
 double yAngleError1 = 0;//length * RandomNG::uniform(-AngleError,AngleError);

 //Jaw 2
 //Random x,y
 double xOffsetError2 = 0;
 double yOffsetError2 = 0;//RandomNG::uniform(-PositionError,PositionError);
 //cout << "xOffsetError1" << "\t" << xOffsetError1 << "\t" << "xOffsetError2" << "\t" << xOffsetError2 << endl;

 //Random theta1, theta2 - small angle approx
 double xAngleError2 = 0;
 double yAngleError2 = 0;//length * RandomNG::uniform(-AngleError,AngleError);
 #endif

 //ENTRANCE
 //We will want to calculate the position of the left and right jaws
 //First
 double gap_x_entrance = collimator_aperture_width_entrance/2;
 double x_rot = (x_orbit * cos(-collimator_aperture_tilt)) - (y_orbit * sin(-collimator_aperture_tilt));
 x_rot = x_orbit;
 cout << "Whole Offset Displacement TDCQA" << "\t" << wholeOffsetDisp << endl;
 double x1 = x_rot + gap_x_entrance + xOffsetError1 ;
 double x2 = x_rot - gap_x_entrance + xOffsetError2  + wholeOffsetDisp;

 //Second
 double gap_y_entrance = collimator_aperture_height_entrance/2;
 double y_rot = (x_orbit * sin(-collimator_aperture_tilt)) + (y_orbit * cos(-collimator_aperture_tilt));
 y_rot = y_orbit;

 double y1 = y_rot + gap_y_entrance + yOffsetError1;
 double y2 = y_rot - gap_y_entrance + yOffsetError2;

 //EXIT
 //First
 //double gap_x_exit = collimator_aperture_width_exit/2;
 double gap_x_exit = collimator_aperture_width_entrance/2;
 double x_rot_exit = (x_orbit_exit * cos(-collimator_aperture_tilt)) - (y_orbit_exit * sin(-collimator_aperture_tilt));
 x_rot_exit = x_orbit_exit;

 //double xx1 = x_rot_exit + gap_x_exit + xOffsetError1 + xAngleError1;
 //double xx2 = x_rot_exit - gap_x_exit + xOffsetError2 + xAngleError2;
 //cout << "xx1" << "\t" << xx1 << "\t" << "xx2" << "\t" << xx2  << endl;

 double xx1 = x_rot_exit + gap_x_exit + xOffsetError1 + xAngleError1  ;
 double xx2 = x_rot_exit - gap_x_exit + xOffsetError2 + xAngleError2  + wholeOffsetDisp;

 //Second
 double gap_y_exit = collimator_aperture_height_exit/2;
 double y_rot_exit = (x_orbit_exit * sin(-collimator_aperture_tilt)) + (y_orbit_exit * cos(-collimator_aperture_tilt));
 y_rot_exit = y_orbit_exit;

 //double yy1 = yy_rot + gap_y_exit;
 //double yy2 = yy_rot - gap_y_exit;

 //double yy1 = y_rot_exit + gap_y_exit + yOffsetError1 + yAngleError1;
 //double yy2 = y_rot_exit - gap_y_exit + yOffsetError2 + yAngleError2;
 double yy1 = y_rot_exit + gap_y_exit + yOffsetError1;
 double yy2 = y_rot_exit - gap_y_exit + yOffsetError2;

 //cout << (CMapit->second)->GetName() << "\t" << x_rot << "\t" << x_rot_exit << "\t" << gap_x_entrance << "\t" << gap_x_exit << "\t" << collimator_aperture_width_entrance << "\t" <<collimator_aperture_width_exit << endl;
 cout << (CMapit->second)->GetName() << "\t\t" << "x1" << "\t\t\t" << "x2" << "\t\t\t" << "xx1" << "\t\t\t" << "xx2" << "\t\t" << "gap_x_entrance" << "\t\t" << "gap_y_entrance" << endl;
 cout << (CMapit->second)->GetName() << "\t" <<  x1  << "\t" <<  x2  << "\t" <<  xx1  << "\t" <<  xx2  << "\t" <<  gap_x_entrance  << "\t" <<  gap_y_entrance << endl;
 //cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<gap_x_entrance << "\t" <<gap_y_entrance << endl;
 //cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<(x1-xx1)*2*1000*1000 << endl;
 //cout << (CMapit->second)->GetName() << "\t" << x1*2 <<"\t" <<  xx1*2 << "\t" <<x2*2 << "\t" << xx2*2<< endl;

 //                      double xj1 = max(x1,xx1);   //+ve values
 //                      double xj2 = min(x2,xx2);   //-ve values

 //                      double yj1 = max(y1,yy1);   //+ve values
 //                      double yj2 = min(y2,yy2);   //-ve values

 double x_size_entrance = (x1 - x2);
 double y_size_entrance = (y1 - y2);

 double x_pos_entrance = (x1 + x2)/2;
 double y_pos_entrance = (y1 + y2)/2;

 double x_size_exit = (xx1 - xx2);
 double y_size_exit = (yy1 - yy2);

 double x_pos_exit = (xx1 + xx2)/2;
 double y_pos_exit = (yy1 + yy2)/2;
 cout << (CMapit->second)->GetName() << "\t\t" << "x_size_entrance" << "\t\t" << "x_pos_entrance" << "\t\t" << "x_size_exit" << "\t\t" << "x_pos_exit" << endl;
 cout << (CMapit->second)->GetName() << "\t" <<  x_size_entrance  << "\t" <<  x_pos_entrance  << "\t" <<  x_size_exit  << "\t" <<  x_pos_exit   << endl;

 //x_size += xOffsetError1;
 //y_size += yOffsetError1;

 CollimatorAperture* app=new CollimatorAperture(x_size_entrance,y_size_entrance,\
                        collimator_aperture_tilt,collimator_material,length,x_pos_entrance,y_pos_entrance);

 app->SetExitXOffset(x_pos_exit);    //Horizontal
 app->SetExitYOffset(y_pos_exit);    //Vertical

 app->SetExitWidth(x_size_exit); //Horizontal
 app->SetExitHeight(y_size_exit);    //Vertical

 //Set the aperture for collimation
 (CMapit->second)->SetAperture(app);
 }
 else
 {
 std::cerr << "Unsupported collimator configuration in CollimatorDatabase.cpp" << std::endl;
 exit(1);
 }

 }
 }
 }
 }
 }

 */

/*
 void CollimatorDatabase::EnableOneSideJawTCDQA(bool flg)
 {
 OneSideJawTCDQA =flg;
 }
 */
void CollimatorDatabase::SelectImpactFactor(string pcoll, double impact)
{
	PrimaryCollimator = pcoll;
	RequestedImpactFactor = impact;     //Impact factor in m
}

void CollimatorDatabase::SetLogFile(ostream& os)
{
	log = &os;
}

void CollimatorDatabase::EnableLogging(bool flg)
{
	logFlag = flg;
}

void CollimatorDatabase::SetErrorLogFile(ostream& os)
{
	ErrorLog = &os;
}

void CollimatorDatabase::EnableErrorLogging(bool flg)
{
	ErrorLogFlag = flg;
}

void CollimatorDatabase::MatchBeamEnvelope(bool flg)
{
	EnableMatchBeamEnvelope = flg;
}

void CollimatorDatabase::EnableJawAlignmentErrors(bool flg)
{
	JawAlignmentErrors = flg;
}

void CollimatorDatabase::SetJawPositionError(double Error)
{
	PositionError = Error;
}

void CollimatorDatabase::SetJawAngleError(double Error)
{
	AngleError = Error;
}

void CollimatorDatabase::OutputFlukaDatabase(std::ostream* os)
{
	(*os)
			<< "# ID\tname\tangle[rad]\tbetax[m]\tbetay[m]\thalfgap[m]\tMaterial\tLength[m]\tsigx[m]\tsigy[m]\ttilt1[rad]\ttilt2[rad]\tnsig"
			<< endl;
	for (vector<FlukaData*>::iterator its = StoredFlukaData.begin(); its != StoredFlukaData.end(); ++its)
	{
		(*os) << setw(6) << left << (*its)->id_coll;
		(*os) << setw(20) << left << (*its)->name;
		//~ (*os) << setw(20) << left << (*its)->position;
		(*os) << setw(12) << left << (*its)->angle;
		(*os) << setw(12) << left << (*its)->beta_x;
		(*os) << setw(12) << left << (*its)->beta_y;
		(*os) << setw(12) << left << (*its)->half_gap;
		(*os) << setw(6) << left << (*its)->material;
		(*os) << setw(12) << left << (*its)->length;
		(*os) << setw(20) << left << (*its)->sig_x;
		(*os) << setw(20) << left << (*its)->sig_y;
		(*os) << setw(12) << left << (*its)->j1_tilt;
		(*os) << setw(12) << left << (*its)->j2_tilt;
		(*os) << setw(12) << left << (*its)->n_sig;
		(*os) << endl;
	}
}

