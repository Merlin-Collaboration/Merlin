/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: (c) a full 6D bunch with wakefields / beam loading
 */

#include "merlin_config.h"

#include <iostream>
#include <fstream>

#include "NumericalConstants.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"
#include "RandomNG.h"

#include "SMPTracker.h"
#include "SMPBunchConstructor.h"
#include "SMPWakeFieldProcess.h"
#include "CorrectorDipoles.h"
#include "TComponentFrame.h"
#include "model_construction.h"
#include "TrackingOutput.h"
#include "QuadReferenceOutput.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace SMPTracking;

// linac tracking benchmarking
// ------------------------------------------
// N. Walker 22.02.2006
//
// This example is based on a code benchmarking exercise.
//
// The two exercises performed here are
//
// 1. track a betatron oscillation along the supplied linac
//    (lattice/tesla_linac.xtff) with an initial vertical
//    offset of 5um as
//    (a) a single ray
//    (b) a full 6D bunch with no wakefields / beam loading
//    (c) a full 6D bunch with wakefields / beam loading
//
// 2. Apply vertical misalignment to the linac lattice as specified
//    in the file lattice/nick23p4_misxy_1.txt. Adjust the vertical
//    dipole correctors as specified in lattice/nick23p4_misxy_ycor_1.txt.
//    Track a full 6D bunch with wakefields and beam loading (no initial offset).
//
// This example uses a Sliced MacroParticle (SMP) tracking for the linac simulation.
// Output is via a purpose written SimulationOutput class (TrackingOuput), which
// calculates various beam parameters at each BPM location.
//
// Output files:
//         quad-reference.data
//         single-particle.data
//         full-bunch-no-WF.data
//         full-bunch-with-WF.data
//         exercise2.data
//

///////////////////////////////////////////////////////////
// forward declarations
///////////////////////////////////////////////////////////

// Tracks the sliced macroparticle bunch through the linac, with or
// without wakefields (transvers and longitudinal).
void PerformTracking(const string& fname, AcceleratorModel::Beamline& linac, SMPBunch* bunch, bool inc_wf,
	SimulationOutput* simout = nullptr);

// Inputs the standard misalignment data and corrector settings
// and adjusts the lattice accordingly (exercise 2)
void AdjustLattice(AcceleratorModel& linac);

///////////////////////////////////////////////////////////
// MAIN PROGRAM
///////////////////////////////////////////////////////////

int main()
{
	// Initialise the random number generated FIRST!
	RandomNG::init(1234);

	// Construct model
	string paths[] = {"../ILCML/lattice/tesla_linac.xtff", "ILCML/lattice/tesla_linac.xtff",
					  "MerlinExamples/ILCML/lattice/tesla_linac.xtff"};

	string lattice_path;
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path = paths[i];
			break;
		}
	}
	cout << "Lattice " << lattice_path << endl;
	pair<AcceleratorModel*, BeamData*> mb = ConstructModel(lattice_path);
	AcceleratorModel* model = mb.first;
	BeamData* beam0 = mb.second;

	AcceleratorModel::Beamline linac = model->GetBeamline();
	SMPBunch* bunch;
	BeamData beam;

	// First task is to output the reference energy and quadrupole
	// fields in the lattice.
	// --------------------------------------------------
	cout << "Calculating quadrupole reference energy..." << flush;
	QuadReferenceOutput* quadRefOutput = new QuadReferenceOutput("quad-reference.data", 5.0, 1.40e+13);
	bunch = SMPBunchConstructor(*beam0, 31, 11).ConstructSMPBunch();
	PerformTracking("dummy", linac, bunch, true, quadRefOutput);
	cout << "done" << endl;

	// Benchmarking Exercise 1
	// --------------------------------------------------

	// Set the required initial beam offset for the oscillation
	beam0->y0 = 5.0 * micrometer;

	cout << "tracking single particle equivalent..." << flush;
	beam = *beam0;
	beam.emit_x = beam.emit_y = beam.sig_dp = beam.sig_z = 0;
	bunch = SMPBunchConstructor(beam, 1, 1).ConstructSMPBunch();
	PerformTracking("single-particle.data", linac, bunch, false);
	cout << "done" << endl;

	cout << "tracking 6D bunch, no WF..." << flush;
	beam = *beam0;
	bunch = SMPBunchConstructor(beam, 31, 11).ConstructSMPBunch();
	PerformTracking("full-bunch-no-WF.data", linac, bunch, false);
	cout << "done" << endl;

	cout << "tracking 6D bunch, with WF..." << flush;
	bunch = SMPBunchConstructor(beam, 31, 11).ConstructSMPBunch();
	PerformTracking("full-bunch-with-WF.data", linac, bunch, true);
	cout << "done" << endl;

	// Benchmarking Exercise 2
	// --------------------------------------------------
	AdjustLattice(*model);
	beam.y0 = 0;
	bunch = SMPBunchConstructor(beam, 31, 11).ConstructSMPBunch();
	PerformTracking("exercise2.data", linac, bunch, true);
	cout << "done" << endl;

	// clean up
	delete model;
	delete beam0;

	return 0;
}

void PerformTracking(const string& fname, AcceleratorModel::Beamline& linac, SMPBunch* bunch, bool inc_wf,
	SimulationOutput* simout)
{
	SMPTracking::SMPTracker tracker(linac);

	SimulationOutput* sout = nullptr;
	if(simout)
	{
		sout = simout;
	}
	else
	{
		TrackingOutput* trackOut = new TrackingOutput(fname);
		trackOut->AddIdentifier("BPM.*");
		trackOut->output_final = true;
		trackOut->output_initial = true;
		trackOut->output_all = false;
		sout = trackOut;
	}

	tracker.SetOutput(sout);

	if(inc_wf)
	{
		WakeFieldProcess* wf = new WakeFieldProcess(1);
		wf->ApplyImpulseAt(WakeFieldProcess::atCentre);
		tracker.AddProcess(wf);
	}

	tracker.Track(bunch);

	delete bunch;
	delete sout;
}

void AdjustLattice(AcceleratorModel& linacModel)
{
	string paths[] = {"../ILCML/lattice/nick23p4_misxy_1.txt", "ILCML/lattice/nick23p4_misxy_1.txt",
					  "MerlinExamples/ILCML/lattice/nick23p4_misxy_1.txt"};

	string lattice_path;
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path = paths[i];
			break;
		}
	}
	// Read in misalignment file
	ifstream ifs(lattice_path);
	if(!ifs)
	{
		cerr << "problems openning file lattice/nick23p4_misxy_1.txt" << endl;
		abort();
	}
	string s;
	getline(ifs, s); // ignore first line

//	AcceleratorModel::Beamline linac = linacModel.GetBeamline();

	vector<ComponentFrame*> frames;
	linacModel.ExtractComponents("TWRFStructure.*|Quadrupole.*|BPM.*", frames);
	cout << "\n" << frames.size() << " components found" << endl;

	for(size_t n = 0; n < frames.size(); n++)
	{
		double z, x, y, xa, ya, tilt;
		ifs >> z >> x >> y >> xa >> ya >> tilt;
		ifs >> s;  // class
		ifs >> s;  // name

		if(!ifs)
		{
			cerr << "Error reading file lattice/nick23p4_misxy_1.txt" << endl;
			abort();
		}

		// check name
		if((*frames[n]).GetComponent().GetName() != s)
		{
			cerr << "name mismatch: ";
			cerr << s << " != " << (*frames[n]).GetComponent().GetName() << endl;
			abort();
		}

		// set alignment
		frames[n]->Translate(x, y, 0);
		frames[n]->RotateY(xa);
		frames[n]->RotateX(-ya); // sign?
		frames[n]->RotateZ(tilt);
	}

	// YCOR settings (assumed units are tesla.meter)
	string path2[] = {"../ILCML/lattice/nick23p4_misxy_ycor_1.txt", "ILCML/lattice/nick23p4_misxy_ycor_1.txt",
					  "MerlinExamples/ILCML/lattice/nick23p4_misxy_ycor_1.txt"};

	string lattice_path2;
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path2 = paths[i];
			break;
		}
	}
	ifstream ifs1(lattice_path2);
	if(!ifs1)
	{
		cerr << "problems openning file lattice/nick23p4_misxy_ycor_1.txt" << endl;
		abort();
	}

	vector<TComponentFrame<YCor>*> ycors;
	linacModel.ExtractTypedComponents(ycors);
	cout << "\n" << ycors.size() << " correctors found" << endl;

	for(size_t n = 0; n < ycors.size(); n++)
	{
		double yc, z;
		ifs1 >> z >> yc;
		if(ifs1)
		{
			cerr << "error reading file lattice/nick23p4_misxy_ycor_1.txt" << endl;
			abort();
		}
		// cout<<z<<" "<<ycors[n]->GetPosition()<<endl;

		YCor& ycor = ycors[n]->GetComponent();
		ycor.SetFieldStrength(yc);
	}
}
