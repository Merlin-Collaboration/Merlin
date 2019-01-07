/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <map>
#include <set>
#include <ctime>

#include "ParticleDistributionGenerator.h"
#include "HaloParticleDistributionGenerator.h"
#include "BunchFilter.h"
#include "ParticleTracker.h"
#include "ParticleBunchTypes.h"
#include "MADInterface.h"
#include "RandomNG.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "CollimateProtonProcess.h"
#include "ScatteringProcess.h"
#include "ScatteringModelsMerlin.h"
#include "CollimatorDatabase.h"
#include "MaterialDatabase.h"
#include "ApertureConfiguration.h"
#include "Dispersion.h"
#include "CollimatorAperture.h"
#include "LossMapCollimationOutput.h"
#include "NANCheckProcess.h"

using namespace std;
using namespace PhysicalUnits;

/*
 * Note this test should fail occasionally. Actually it fails quite often.
 * By default it only uses 1000 particles, which is too few. 10k or 100k
 * are really needed for a more reliable test, but this makes it slow.
 *
 * Compute the loss map for the nominal LHC lattice, and compare with
 * a pre-computed version.
 *
 * arguments:
 *    collimation_test nparticles seed
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

enum loss_map_mode_t
{
	HORIZONTAL_LOSS,
	VERTICAL_LOSS

};

bool SortComponent(const AcceleratorComponent* first, const AcceleratorComponent* last)
{
	return first->GetComponentLatticePosition() < last->GetComponentLatticePosition();
}

vector<AcceleratorComponent*> SortAcceleratorModel(AcceleratorModel* model)
{
	vector<AcceleratorComponent*> elements;
	model->ExtractTypedElements(elements, "*");

	//Now sort the elements in the appropriate lattice order
	sort(elements.begin(), elements.end(), SortComponent);
	return elements;
}

int FindElementLatticePosition(string RequestedElement, AcceleratorModel* model)
{
	vector<AcceleratorComponent*> elements = SortAcceleratorModel(model);
	size_t nelm = elements.size();
	for(size_t n = 0; n < nelm; n++)
	{
		if(elements[n]->GetName() == RequestedElement)
		{
			cout << "Found " << RequestedElement << " at " << n << " of " << nelm << endl;
			return n;
		}
	}
	return 0;
}

int main(int argc, char* argv[])
{
	int seed = 0;
	int npart = 10000;
	//~ int npart = 100000;
	int nturns = 20;

	//Loss_Map or Merged Collimation
	bool Loss_Map               = 0;
	bool st_scatter             = 0;

	if(argc >= 2)
	{
		seed = atoi(argv[1]);
	}

	if(seed == 0)
	{
		seed = (int) time(nullptr);
	}

	if(argc >= 3)
	{
		npart = atoi(argv[2]);
	}

	cout << "Random Seed: " << seed << endl;
	RandomNG::init(seed);

	// find directories
	string log_dir = "";
	string input_data_dir = "";
	string result_dir = "";

	string paths[] = {"../", "", "MerlinTests/"};
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open((paths[i] + "data/collimator.7.0.sigma").c_str());
		if(test_file)
		{
			test_file.close();
			input_data_dir = paths[i] + "data/";
			log_dir = paths[i] + "outputs/";
			result_dir = paths[i] + "outputs/";
			break;
		}
	}
	if(log_dir == "")
	{
		cout << "Could not find data directory. Try running from cmake build directory, or the executable directory"
			 << endl;
		exit(1);
	}

	double beam_energy = 7000.0;

	cout << "npart=" << npart << " nturns=" << nturns << endl;

	loss_map_mode_t loss_map_mode = HORIZONTAL_LOSS;
//	loss_map_mode_t loss_map_mode = VERTICAL_LOSS;

	string start_element;
	switch(loss_map_mode)
	{
	case HORIZONTAL_LOSS:
		start_element = "TCP.C6L7.B1";  //HORIZONTAL COLLIMATOR (x)
		break;
	case VERTICAL_LOSS:
		start_element = "TCP.D6L7.B1";  //VERTICAL COLLIMATOR (y)
		break;
	}

	double beamcharge = 1.1e11;
	double normalized_emittance = 3.5e-6;
	double gamma = beam_energy / PhysicalConstants::ProtonMassMeV / PhysicalUnits::MeV;
	double beta = sqrt(1.0 - (1.0 / pow(gamma, 2)));
	double emittance = normalized_emittance / (gamma * beta);

	MaterialDatabase* mat = new MaterialDatabase();
	CollimatorDatabase* collimator_db = new CollimatorDatabase(input_data_dir + "collimator.7.0.sigma", mat, true);

	//	ACCELERATOR MODEL LOADING
	MADInterface* myMADinterface;
	myMADinterface = new MADInterface(input_data_dir + "twiss.7.0tev.b1_new.tfs", beam_energy);
	myMADinterface->TreatTypeAsDrift("RFCAVITY");

	//Build accelerator model
	AcceleratorModel* model = myMADinterface->ConstructModel();

	LatticeFunctionTable* twiss = new LatticeFunctionTable(model, beam_energy);
	twiss->AddFunction(1, 6, 3);
	twiss->AddFunction(2, 6, 3);
	twiss->AddFunction(3, 6, 3);
	twiss->AddFunction(4, 6, 3);
	twiss->AddFunction(6, 6, 3);

	// find twiss
	double bscale1 = 1e-22;
	while(true)
	{
		twiss->ScaleBendPathLength(bscale1);
		twiss->Calculate();
		if(!std::isnan(twiss->Value(1, 1, 1, 0)))
		{
			break;
		}
		bscale1 *= 2;
	}

	// FLAG to set an automatically matching between beam envelope and collimator taper
	collimator_db->MatchBeamEnvelope(false);
	collimator_db->EnableJawAlignmentErrors(false);

	//Collimator rms error on gap size 0.1 sigma, jaw angle error with respect the beam envelope 200 microradiant
	collimator_db->SetJawPositionError(0.0 * nanometer); // This is actually the variance of the error alignment of 0.1 sigma i.e. 35.3 micronmeter.... it should be in nm^2
	collimator_db->SetJawAngleError(0.0 * microradian); // rms error 200 microradiant
	collimator_db->SelectImpactFactor(start_element, 1.0e-6);

	double impact;
	//Setup the collimator jaws to appropriate sizes and
	try
	{
		impact = collimator_db->ConfigureCollimators(model, emittance, emittance, twiss);
	}
	catch(exception& e)
	{
		std::cout << "Exception caught: " << e.what() << std::endl;
		exit(1);
	}
	if(std::isnan(impact))
	{
		cerr << "Impact is nan" << endl;
		exit(1);
	}
	cout << "Impact factor number of sigmas: " << impact << endl;
	delete collimator_db;

	ApertureConfiguration* apc = new ApertureConfiguration(input_data_dir + "Aperture_B1_6p5TeV_2016.tfs");

	apc->ConfigureElementApertures(model);
	cout << "aperture load finished" << endl << "start twiss" << endl;
	delete apc;

	//Calculate Dispersion
	Dispersion* disp = new Dispersion(model, beam_energy);
	int start_element_number = FindElementLatticePosition(start_element.c_str(), model);
	disp->FindDispersion(start_element_number);

	//      BEAM SETTINGS
	BeamData mybeam;

	//Default values are 0.0
	//The charge of the particles in the beam.
	//  <0 for electrons, >0 for positrons/protons.
	mybeam.charge = beamcharge / npart;
	mybeam.p0 = beam_energy;
	mybeam.beta_x = twiss->Value(1, 1, 1, start_element_number) * meter;
	mybeam.beta_y = twiss->Value(3, 3, 2, start_element_number) * meter;
	mybeam.alpha_x = -twiss->Value(1, 2, 1, start_element_number);
	mybeam.alpha_y = -twiss->Value(3, 4, 2, start_element_number);

	//Dispersion
	mybeam.Dx = disp->Dx;
	mybeam.Dy = disp->Dy;
	//mybeam.Dy=0;
	mybeam.Dxp = disp->Dxp;
	mybeam.Dyp = disp->Dyp;

	mybeam.emit_x = impact * impact * emittance * meter;
	impact = 1;
	mybeam.emit_y = impact * impact * emittance * meter;
	mybeam.sig_z = 0.0;

	//Beam centroid
	mybeam.x0 = twiss->Value(1, 0, 0, start_element_number);
	mybeam.xp0 = twiss->Value(2, 0, 0, start_element_number);
	mybeam.y0 = twiss->Value(3, 0, 0, start_element_number);
	mybeam.yp0 = twiss->Value(4, 0, 0, start_element_number);
	mybeam.ct0 = twiss->Value(5, 0, 0, start_element_number);

	mybeam.sig_dp = 0.0;

	//X-Y coupling
	mybeam.c_xy = 0.0;
	mybeam.c_xyp = 0.0;
	mybeam.c_xpy = 0.0;
	mybeam.c_xpyp = 0.0;

	delete disp;

	//	BUNCH SETTINGS
	int node_particles = npart;

	vector<Collimator*> TCP;
	int siz = model->ExtractTypedElements(TCP, start_element);
	Aperture *ap = (TCP[0])->GetAperture();
	if(!ap)
	{
		cout << "Could not get tcp ap" << endl;
		abort();
	}
	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
	if(!CollimatorJaw)
	{
		cout << "Could not cast" << endl;
		abort();
	}

	double h_offset = twiss->Value(1, 0, 0, start_element_number);
	double JawPosition = CollimatorJaw->GetFullEntranceWidth() / 2.0;
	HorizontalHaloParticleBunchFilter *hFilter = new HorizontalHaloParticleBunchFilter();
	hFilter->SetHorizontalLimit(JawPosition);
	hFilter->SetHorizontalOrbit(h_offset);

	// dist 1 is halo in 1 plane, zero in other
	// dist 2 is halo in 1 plane, gauss in other
	ProtonBunch* myBunch;
	switch(loss_map_mode)
	{
	case HORIZONTAL_LOSS:
		myBunch = new ProtonBunch(node_particles, HorizonalHalo2ParticleDistributionGenerator(), mybeam, hFilter);
		break;
	case VERTICAL_LOSS:
		myBunch = new ProtonBunch(node_particles, VerticalHalo2ParticleDistributionGenerator(), mybeam, hFilter);
		break;
	}

	myBunch->SetMacroParticleCharge(mybeam.charge);

	if(Loss_Map && st_scatter)
	{
		myBunch->EnableScatteringPhysics(ProtonBunch::SixTrack);
	}
	else if(Loss_Map && !st_scatter)
	{
		myBunch->EnableScatteringPhysics(ProtonBunch::Merlin);
	}

	//	PARTICLE TRACKER
	AcceleratorModel::RingIterator bline = model->GetRing(start_element_number);
	ParticleTracker* tracker = new ParticleTracker(bline, myBunch);
	tracker->SetLogStream(std::cout);

	auto nancheck = new NANCheckProcess;
	nancheck->SetDetailed(0);
	nancheck->SetHaltNAN(1);
	tracker->AddProcess(nancheck);

	//	COLLIMATION SETTINGS
	std::ostringstream filename;
	filename << "lhc_collimation_test_lossmap_" << npart << ".dat";

	LossMapCollimationOutput* myLossOutput = new LossMapCollimationOutput(tencm);
	ScatteringModel* myScatter;
	if(Loss_Map)
	{
		CollimateParticleProcess* myCollimateProcess;

		myCollimateProcess = new CollimateParticleProcess(2, 4);
		myCollimateProcess->ScatterAtCollimator(true);

		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);
		myCollimateProcess->SetCollimationOutput(myLossOutput);
		tracker->AddProcess(myCollimateProcess);
	}
	else
	{
		CollimateProtonProcess* myCollimateProcess;

		myCollimateProcess = new CollimateProtonProcess(2, 4);
		myCollimateProcess->ScatterAtCollimator(true);

		if(st_scatter)
		{
			myScatter = new ScatteringModelSixTrack;
		}
		else
		{
			myScatter = new ScatteringModelMerlin;
		}

		myCollimateProcess->SetScatteringModel(myScatter);

		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);
		myCollimateProcess->SetCollimationOutput(myLossOutput);
		tracker->AddProcess(myCollimateProcess);
	}

	//	TRACKING RUN
	for(int turn = 1; turn <= nturns; turn++)
	{
		cout << "Turn " << turn << "\tParticle number: " << myBunch->size() << endl;
		tracker->Track(myBunch);
		if(myBunch->size() <= 1)
		{
			break;
		}
	}

	myLossOutput->Finalise();

	ofstream* col_output = new ofstream(filename.str().c_str());
	if(!col_output->good())
	{
		std::cerr << "Could not open collimation loss file" << std::endl;
		exit(EXIT_FAILURE);
	}
	*col_output << "# From merlin npart=" << npart << " nturns=" << nturns << endl;
	myLossOutput->Output(col_output);

	col_output->flush();
	col_output->close();
	cout << "Wrote: " << filename.str() << endl;
	delete col_output;

	cout << "npart: " << npart << endl;
	cout << "left: " << myBunch->size() << endl;
	cout << "absorbed: " << npart - myBunch->size() << endl;

	delete mat;
	delete myBunch;
	delete twiss;
	delete model;
	delete myMADinterface;
}
