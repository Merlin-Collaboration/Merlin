/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <unistd.h>
#include <cstdint>
#include <string>

#include "Components.h"
#include "CollimatorAperture.h"
#include "AcceleratorModelConstructor.h"

#include "ParticleTracker.h"
#include "ParticleBunchTypes.h"

#include "CollimateParticleProcess.h"
#include "CollimateProtonProcess.h"
#include "ScatteringModelsMerlin.h"
#include "MaterialData.h"
#include "NANCheckProcess.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

#include "RandomNG.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

/*
 * Note this test should fail occasionally
 *
 * see cu50_test.py for details
 *
 * The time (seconds since epoch) is used as the seed unless an int is
 * passed as the first argument.
 *
 *
 */

int main(int argc, char* argv[])
{
	int seed = 0;
	if(argc >= 2)
	{
		seed = atoi(argv[1]);
	}

	if(seed == 0)
	{
		seed = (int) time(nullptr);
	}

	//Number of particles
	//const int npart = 1e9;
	uint64_t npart = 1e6;
	if(argc >= 3)
	{
		npart = std::stoull(argv[2]);
	}

	bool scatter_mode_sixtrack = 0;

	if(argc >= 4)
	{
		for(int i = 3; i < argc; i++)
		{
			cout << "opt " << argv[i] << endl;
			if(strcmp(argv[i], "sixtrack") == 0)
			{
				scatter_mode_sixtrack = 1;
			}
		}
	}

	cout << "Seed: " << seed << endl;
	RandomNG::init(seed);
	/*********************************************************************
	 *	GENERAL SETTINGS
	 *********************************************************************/
	bool output_final_bunch     = 0;

	//Beam energy (GeV) 7000,3500,450 etc
	//double beam_energy = 7000.0;
	const double beam_energy = 7000.0;

	uint64_t particles_left = npart;
	const uint64_t max_particles_per_bunch = 1e6; // batch, to avoid high mem usage

	const size_t nbins = 100;

	double bin_mins[5], bin_maxs[5];
	bin_mins[0] = -50e-6;
	bin_maxs[0] = 50e-6;   // x
	bin_mins[1] = -100e-6;
	bin_maxs[1] = 100e-6; //xp
	bin_mins[2] = -50e-6;
	bin_maxs[2] = 50e-6;   //y
	bin_mins[3] = -100e-6;
	bin_maxs[3] = 100e-6; //yp
	bin_mins[4] = 1e-6;
	bin_maxs[4] = 1e-1;      //dp

	uint64_t hists[5][nbins + 2] = {{0}};

	/*********************************************************************
	 *	ACCELERATOR MODEL LOADING
	 *********************************************************************/

	StandardMaterialData* mat = new StandardMaterialData();
	if(scatter_mode_sixtrack)
		mat->UseSixTrackValues();
	cout << mat << endl;
	MaterialProperties* CollimatorMaterial = mat->property["Cu"];

	AcceleratorModelConstructor* construct = new AcceleratorModelConstructor();
	double length = 0.5;
	Collimator* TestCol = new Collimator("TestCollimator", length);
	TestCol->SetMaterialProperties(CollimatorMaterial);

	CollimatorAperture* app = new CollimatorAperture(2, 2, 0, length, 0, 0);
	app->SetExitWidth(app->GetFullEntranceWidth());      //Horizontal
	app->SetExitHeight(app->GetFullEntranceHeight());    //Vertical

	//Set the aperture for collimation
	TestCol->SetAperture(app);

	construct->AppendComponent(*TestCol);

	AcceleratorModel* model = construct->GetModel();

	ostringstream bunch_output_file2;
	ofstream bunch_output2;
	if(output_final_bunch)
	{
		bunch_output_file2 << "Bunch/LM_M_final.txt";
		bunch_output2.open(bunch_output_file2.str().c_str());
	}

	/*********************************************************************
	 *	PARTICLE TRACKER
	 *********************************************************************/
	ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1);
	AcceleratorModel::RingIterator bline = model->GetRing();
	ParticleTracker* tracker = new ParticleTracker(bline, myBunch);

	auto nancheck = new NANCheckProcess;
	nancheck->SetDetailed(0);
	nancheck->SetHaltNAN(1);
	tracker->AddProcess(nancheck);

	/*********************************************************************
	 *	COLLIMATION SETTINGS
	 *********************************************************************/
	ScatteringModel* myScatter = new ScatteringModel();

	CollimateProtonProcess* myCollimateProcess = new CollimateProtonProcess(2, 4);
	if(scatter_mode_sixtrack)
	{
// RJB		myScatter = new ScatteringModelSixTrack;
		myScatter->Processes[1] = new SixTrackRutherford();
		myScatter->Processes[2] = new SixTrackElasticpn();
		myScatter->Processes[3] = new SixTrackSingleDiffractive();
	}
//	else
//	{
//		myScatter = new ScatteringModelMerlin;
//	}
	myCollimateProcess->SetScatteringModel(myScatter);
	stringstream loststr;

	myCollimateProcess->ScatterAtCollimator(true);

	// Sets maximum allowed loss percentage at a single collimator.
	myCollimateProcess->SetLossThreshold(101.0);

	//sets process log stream, nullptr to disable. aka, what col_output is above
	myCollimateProcess->SetLogStream(nullptr);

	//Add Collimation process to the tracker.
	myCollimateProcess->SetOutputBinSize(length);
	tracker->AddProcess(myCollimateProcess);

	double y_offset = 1.0 + 1e-6;
	while(particles_left > 0)
	{
		cout << "Particles_left " << particles_left << " of " << npart << endl;
		/*********************************************************************
		 *      BEAM SETTINGS
		 *********************************************************************/
		size_t this_npart = particles_left;
		if(this_npart > max_particles_per_bunch)
		{
			//cout << "Run subset of " << max_particles_per_bunch << endl;
			this_npart = max_particles_per_bunch;
		}

		particles_left -= this_npart;
		//cout << "this_npart " << this_npart << " new particles_left " << particles_left << endl;
		myBunch->clear();
		Particle p(0);

		p.y() = y_offset;
		for(size_t i = 0; i < this_npart; i++)
		{
			myBunch->AddParticle(p);
		}

		/*********************************************************************
		 *	Tracking
		 *********************************************************************/

		cout << "Tracking" << endl;
		tracker->Track(myBunch);
		cout << "Finished.\tParticle number: " << myBunch->size() << endl;
		cout << "npart: " << this_npart << endl;
		cout << "left: " << myBunch->size() << endl;
		cout << "absorbed: " << this_npart - myBunch->size() << endl;

		/*********************************************************************
		 *	Output Final Bunch
		 *********************************************************************/

		if(output_final_bunch)
		{
			myBunch->Output(bunch_output2);
		}
		// Histogramming
		for(PSvectorArray::iterator ip = myBunch->begin(); ip != myBunch->end(); ip++)
		{
			double coords[] = {ip->x(), ip->xp(), ip->y() - y_offset, ip->yp(), -ip->dp()};
			for(int i = 0; i < 5; i++)
			{
				// beware, this can rollover when x is big
				int bin_x = ((coords[i] - bin_mins[i]) / (bin_maxs[i] - bin_mins[i]) * (nbins)) + 1; // +1 because bin zero for outliers
				// so handle end bins, by check against x, not bin
				if(coords[i] < bin_mins[i])
				{
					bin_x = 0;
				}
				if(coords[i] > bin_maxs[i])
				{
					bin_x = nbins + 1;
				}
				hists[i][bin_x] += 1;
			}
		}
	}

	delete myScatter;
	delete myBunch;
	delete tracker;
	delete construct;
	delete mat;
	delete TestCol;
	/*********************************************************************
	 *	Output Final Hist
	 *********************************************************************/
	std::ostringstream filename;
	filename << "cu50_test_hist_" << (scatter_mode_sixtrack ? "st" : "m") << "_" << npart << ".dat";

	std::ofstream out2;
	out2.open(filename.str());
	out2 << "# From merlin npart=" << npart << endl;
	out2 << "# bin x xp y yp dp" << endl;
	for(size_t i = 0; i < nbins + 2; i++)
	{
		out2 << i << " ";
		for(int j = 0; j < 5; j++)
		{
			out2 << hists[j][i] << " ";
		}
		out2 << endl;
	}
	cout << "Wrote " << filename.str() << endl;

	return 0;
}
