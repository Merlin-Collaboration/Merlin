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

#include "Components.h"
#include "CollimatorAperture.h"
#include "AcceleratorModelConstructor.h"

#include "ParticleTracker.h"
#include "ParticleBunchTypes.h"

#include "CollimateParticleProcess.h"
#include "CollimateProtonProcess.h"
#include "MaterialDatabase.h"

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
 * Sends a beam of particles into a 50cm copper block, as described
 * in "TOOLS FOR PREDICTING CLEANING EFFICIENCY IN THE LHC" by
 * R. Assmann et al, PAC2003
 *
 * This is use Merlin's emulation of the SixTrack scattering model.
 *
 * Results are compared to a distribution generated with merlin in
 * 2015, which is roughly in agreement with those by Assmann. Hence
 * this currently checks that the physics have not changed, rather
 * than checking that the physics is correct.
 *
 * The number of bins deviating by more than N sigma are recorded
 * and if the are significantly more that expected (by a normal
 * distribution), the test fails.
 *
 * If the test fails run it a few times. If it consistently fails then
 * there is an actual issue. Bins out by more than 3 sigma are
 * displayed, repeated failure of the same bin would be suspicious.
 *
 * The time (seconds since epoch) is used as the seed unless an int is
 * passed as the first argument.
 *
 *
 */

int main(int argc, char* argv[])
{
	int seed = (int) time(NULL);
	if(argc >= 2)
	{
		seed = atoi(argv[1]);
	}

	cout << "Seed: " << seed << endl;
	RandomNG::init(seed);
	/*********************************************************************
	 *	GENERAL SETTINGS
	 *********************************************************************/
	//Loss_Map or Merged Collimation
	bool Loss_Map               = 0;
	if(Loss_Map)
	{
		std::cout << "LOSSMAP Collimation (ProtonBunch)" << std::endl;
	}
	else
	{
		std::cout << "MERGED Collimation (based on HEL code)" << std::endl;
	}
	bool output_initial_bunch   = 0;
	bool output_final_bunch     = 1;

	//Beam energy (GeV) 7000,3500,450 etc
	//double beam_energy = 7000.0;
	const double beam_energy = 7000.0;

	//Number of particles
	const int npart = 1e6;
	//~ const int npart = 10;

	const size_t nbins = 100;
	const double bin_min_x = -50e-6, bin_max_x = 50e-6;
	const double bin_min_xp = -100e-6, bin_max_xp = 100e-6;
	const double bin_min_dp = 1e-3, bin_max_dp = 1e-1;

	int hist_x[nbins + 2] = {0};
	int hist_xp[nbins + 2] = {0};
	int hist_dp[nbins + 2] = {0};

	// load data file
	int c_hist_x[nbins + 2] = {0};
	int c_hist_xp[nbins + 2] = {0};
	int c_hist_dp[nbins + 2] = {0};

	string paths[] = {"../data/cu50_merlin_sixtrack_hist.dat", "data/cu50_merlin_sixtrack_hist.dat", "MerlinTests/data/cu50_merlin_sixtrack_hist.dat"};
	ifstream dist_file;
	for(size_t i = 0; i < 3; i++)
	{
		dist_file.open(paths[i].c_str());
		if(dist_file)
		{
			break;
		}
	}
	if(!dist_file)
	{
		cout << "Failed to open dist file" << endl;
		exit(1);
	}

	size_t dist_bin;
	double dist_val_x, dist_val_xp, dist_val_dp;
	string line;
	while(getline(dist_file, line))
	{
		if(line == "")
		{
			continue;
		}
		if(line[0] == '#')
		{
			continue;
		}
		istringstream liness(line);
		liness >> dist_bin >> dist_val_x >> dist_val_xp >> dist_val_dp;
		if(dist_bin >= (nbins + 2))
		{
			cout << "Invalid bin (" << dist_bin << ") in cu50_merlin_hist.dat";
			exit(1);
		}
		c_hist_x[dist_bin] = dist_val_x * npart;
		c_hist_xp[dist_bin] = dist_val_xp * npart;
		c_hist_dp[dist_bin] = dist_val_dp * npart;
	}

	/*********************************************************************
	 *	ACCELERATOR MODEL LOADING
	 *********************************************************************/

	MaterialDatabase* mat = new MaterialDatabase();
	Material* CollimatorMaterial = mat->FindMaterial("Cu");

	AcceleratorModelConstructor* construct = new AcceleratorModelConstructor();
	construct->NewModel();
	double length = 0.5;
	Collimator* TestCol = new Collimator("TestCollimator", length);
	TestCol->SetMaterial(CollimatorMaterial);

	CollimatorAperture* app = new CollimatorAperture(2, 2, 0, CollimatorMaterial, length, 0, 0);
	app->SetExitWidth(app->GetFullEntranceWidth());      //Horizontal
	app->SetExitHeight(app->GetFullEntranceHeight());    //Vertical

	//Set the aperture for collimation
	TestCol->SetAperture(app);

	construct->AppendComponent(*TestCol);
	AcceleratorModel* model = construct->GetModel();

	/*********************************************************************
	 *      BEAM SETTINGS
	 *********************************************************************/
	ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1);
	Particle p(0);
	p.y() = 1.0 + 1e-6;
	for(int i = 0; i < npart; i++)
	{
		myBunch->AddParticle(p);
	}

	/*********************************************************************
	 *	Output Initial Bunch
	 *********************************************************************/
	if(output_initial_bunch)
	{
		ostringstream bunch_output_file;
		if(Loss_Map)
		{
			bunch_output_file << "Bunch/LM_ST_initial.txt";
		}
		else
		{
			bunch_output_file << "Bunch/HEL_ST_initial.txt";
		}

		ofstream* bunch_output = new ofstream(bunch_output_file.str().c_str());
		myBunch->Output(*bunch_output);
		delete bunch_output;
	}

	/*********************************************************************
	 *	PARTICLE TRACKER
	 *********************************************************************/
	AcceleratorModel::RingIterator bline = model->GetRing();
	ParticleTracker* tracker = new ParticleTracker(bline, myBunch);

	/*********************************************************************
	 *	COLLIMATION SETTINGS
	 *********************************************************************/
	if(Loss_Map)
	{
		CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(2, 4);
		//myBunch->EnableScatteringPhysics(ProtonBunch::Merlin);
		myBunch->EnableScatteringPhysics(ProtonBunch::SixTrack);
		stringstream loststr;

		myCollimateProcess->ScatterAtCollimator(true);

		// Sets maximum allowed loss percentage at a single collimator.
		myCollimateProcess->SetLossThreshold(101.0);

		//sets process log stream, NULL to disable. aka, what col_output is above
		myCollimateProcess->SetLogStream(NULL);

		//Add Collimation process to the tracker.
		myCollimateProcess->SetOutputBinSize(length);
		tracker->AddProcess(myCollimateProcess);
	}
	else
	{
		CollimateProtonProcess* myCollimateProcess = new CollimateProtonProcess(2, 4);
		ScatteringModel* myScatter = new ScatteringModel;
		myScatter->SetScatterType(0);
		myCollimateProcess->SetScatteringModel(myScatter);
		stringstream loststr;

		myCollimateProcess->ScatterAtCollimator(true);

		// Sets maximum allowed loss percentage at a single collimator.
		myCollimateProcess->SetLossThreshold(101.0);

		//sets process log stream, NULL to disable. aka, what col_output is above
		myCollimateProcess->SetLogStream(NULL);

		//Add Collimation process to the tracker.
		myCollimateProcess->SetOutputBinSize(length);
		tracker->AddProcess(myCollimateProcess);
	}
	/*********************************************************************
	 *	Tracking
	 *********************************************************************/

	cout << "Tracking" << endl;
	tracker->Track(myBunch);
	cout << "Finished.\tParticle number: " << myBunch->size() << endl;
	cout << "npart: " << npart << endl;
	cout << "left: " << myBunch->size() << endl;
	cout << "absorbed: " << npart - myBunch->size() << endl;

	if(0)
	{
		ostringstream bunch_output_file_out;
		bunch_output_file_out << "cu50_test_bunch_out.txt";
		ofstream* bunch_output_out = new ofstream(bunch_output_file_out.str().c_str());
		*bunch_output_out << "#T0 P0 x xp y yp ct dp" << endl;
		myBunch->Output(*bunch_output_out);
	}

	/*********************************************************************
	 *	Output Final Bunch
	 *********************************************************************/
	if(output_final_bunch)
	{
		ostringstream bunch_output_file2;
		if(Loss_Map)
		{
			bunch_output_file2 << "Bunch/LM_ST_final.txt";
		}
		else
		{
			bunch_output_file2 << "Bunch/HEL_ST_final.txt";
		}

		ofstream* bunch_output2 = new ofstream(bunch_output_file2.str().c_str());
		myBunch->Output(*bunch_output2);
		delete bunch_output2;
	}

	// Histogramming
	for(PSvectorArray::iterator ip = myBunch->begin(); ip != myBunch->end(); ip++)
	{

		int bin_x = ((ip->x() - bin_min_x) / (bin_max_x - bin_min_x) * (nbins)) + 1; // +1 because bin zero for outliers
		//cout << "\nbin_x = " << bin_x << " x = " << ip->x() << endl;
		// so handle end bins, by check against x, not bin
		if(ip->x() < bin_min_x)
		{
			bin_x = 0;
		}
		if(ip->x() > bin_max_x)
		{
			bin_x = nbins + 1;
		}
		hist_x[bin_x] += 1;

		int bin_xp = ((ip->xp() - bin_min_xp) / (bin_max_xp - bin_min_xp) * (nbins)) + 1; // +1 because bin zero for outliers
		if(ip->xp() < bin_min_xp)
		{
			bin_xp = 0;
		}
		if(ip->xp() > bin_max_xp)
		{
			bin_xp = nbins + 1;
		}
		hist_xp[bin_xp] += 1;

		int bin_dp = ((-ip->dp() - bin_min_dp) / (bin_max_dp - bin_min_dp) * (nbins)) + 1; // +1 because bin zero for outliers
		if(-ip->dp() < bin_min_dp)
		{
			bin_dp = 0;
		}
		if(-ip->dp() > bin_max_dp)
		{
			bin_dp = nbins + 1;
		}
		hist_dp[bin_dp] += 1;
	}

	delete myBunch;

	//~ cout << "i x xp dp" << endl;
	//~ for (size_t i=0; i<nbins+2; i++){
	//~ cout << i << " " << (double)hist_x[i]/npart << " " << (double)hist_xp[i]/npart <<" " << (double)hist_dp[i]/npart << endl;
	//~ }

	/*********************************************************************
	 *	Output Final Hist
	 *********************************************************************/
	if(Loss_Map)
	{
		std::ofstream out2("Bunch/LM_S_hist.txt", std::ofstream::out);
		for(size_t i = 0; i < nbins + 2; i++)
		{
			out2 << i << " " << (double) hist_x[i] / npart << " " << (double) hist_xp[i] / npart << " " << (double) hist_dp[i] / npart << endl;
		}
	}
	else
	{
		std::ofstream out2("Bunch/HEL_S_hist.txt", std::ofstream::out);
		for(size_t i = 0; i < nbins + 2; i++)
		{
			out2 << i << " " << (double) hist_x[i] / npart << " " << (double) hist_xp[i] / npart << " " << (double) hist_dp[i] / npart << endl;
		}
	}

	// Calculate sigma from sqrt(<n>)
	// Count number of bins out by more than n-sigma
	const int ns = 4;
	// https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule, could use erf() in c++11
	const double tails[] = {0, 0.682689492137086, 0.954499736103642, 0.997300203936740,
							0.999936657516334, 0.999999426696856, 0.999999998026825};
	int over_ns[ns + 1] = {0};
	const double error = 0.05; // accept a 5% error in addtion to stat error

	for(size_t i = 0; i < nbins + 2; i++)
	{
		// should use poisson stats or ks test, otherwise bins with small counts can give bad results
		// Currently assume gaussian, don't allow sigma to be below 1.0
		double exp_err_x = max(sqrt(c_hist_x[i]), 1.0);
		double diff_x = hist_x[i] - c_hist_x[i];
		exp_err_x = sqrt(pow(exp_err_x, 2) + pow(c_hist_x[i] * error, 2));

		double exp_err_xp = max(sqrt(c_hist_xp[i]), 1.0);
		double diff_xp = hist_xp[i] - c_hist_xp[i];
		exp_err_xp = sqrt(pow(exp_err_xp, 2) + pow(c_hist_xp[i] * error, 2));

		double exp_err_dp = max(sqrt(c_hist_dp[i]), 1.0);
		double diff_dp = hist_dp[i] - c_hist_dp[i];
		exp_err_dp = sqrt(pow(exp_err_dp, 2) + pow(c_hist_dp[i] * error, 2));

		for(int j = 1; j < ns + 1; j++)
		{
			if(fabs(diff_x) / exp_err_x > j)
			{
				over_ns[j]++;
				// display the worst bins
				if(j >= 3)
				{
					cout << endl;
					cout << "big error in x # bin, merlin, expected, error/sigma" << endl;
					cout << i << " " << hist_x[i] << " " << c_hist_x[i] << "+-" << exp_err_x << " " << diff_x / exp_err_x << endl;
				}
			}
			if(fabs(diff_xp) / exp_err_xp > j)
			{
				over_ns[j]++;
				// display the worst bins
				if(j >= 3)
				{
					cout << endl;
					cout << "big error in xp # bin, merlin, expected, error/sigma" << endl;
					cout << i << " " << hist_xp[i] << " " << c_hist_xp[i] << "+-" << exp_err_xp << " " << diff_xp / exp_err_xp << endl;
				}
			}
			if(fabs(diff_dp) / exp_err_dp > j)
			{
				over_ns[j]++;
				// display the worst bins
				if(j >= 3)
				{
					cout << endl;
					cout << "big error in dp # bin, merlin, expected, error/sigma" << endl;
					cout << i << " " << hist_dp[i] << " " << c_hist_dp[i] << "+-" << exp_err_dp << " " << diff_dp / exp_err_dp << endl;
				}
			}
		}
	}

	// put limit on fraction of bins out by more than n-sigma
	// bad if 30% more than expected
	// allow some leeway for low count bins, by using over_ns[j]-1)/(nbins+2)

	bool bad = 0;
	double tol = 1.3;
	cout << endl << "Bins out by more than";
	for(int j = 1; j < ns + 1; j++)
	{
		cout << endl << "    " << j << " sigma: " << over_ns[j] << " (" << (double) over_ns[j] / (nbins + 2) / 3 * 100  << " %)";
		if((double) over_ns[j] / (nbins + 2) / 3 > (1 - tails[j]) * tol)
		{
			cout << "** bad **";
			bad = 1;
		}
	}
	cout << endl;

	return bad;
}
