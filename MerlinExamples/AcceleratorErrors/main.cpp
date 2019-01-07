/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

// an example how to use the AcceleratorError class
#include "merlin_config.h"
#include "RandomNG.h"
#include "AcceleratorModel.h"
#include "XTFFInterface.h"
#include "BeamData.h"
#include "SMPComponentTracker.h"
#include "SMPTracker.h"
#include "SMPBunch.h"
#include "SMPBunchConstructor.h"

#include "TWRFStructure.h"
#include "AcceleratorSupport.h"
#include "SupportStructure.h"

#include "NumericalConstants.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"

//#include "SplitMagnets.h"
using namespace std;

#include "AcceleratorErrors.h"

#include <iostream>

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace SMPTracking;

pair<AcceleratorModel*, BeamData*> ConstructModel(const string& fname);
void Print(string tag, SMPBunch* bunch);

int main()
{

	// Initialise the random number generator
	RandomNG::init();

	// Construct model
	string paths[] = {"../lattices/ilc_linac_15_250.xtff", "lattices/ilc_linac_15_250.xtff",
					  "MerlinExamples/lattices/ilc_linac_15_250.xtff"};

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

	AcceleratorModel*          model = mb.first;
	BeamData*                  beam  = mb.second;
	AcceleratorModel::Beamline bline = model->GetBeamline();

	// a sub beam line
	AcceleratorModel::Beamline blineFirst = model->GetBeamline(1, 1000);

	// Construct 2 bunches according to beam parameter
	SMPBunch* bunch1 = SMPBunchConstructor(*beam, 31, 11).ConstructSMPBunch();
	SMPBunch* bunch2 = SMPBunchConstructor(*beam, 31, 11).ConstructSMPBunch();

	Print("Initial ", bunch1);

	// a tracker
	SMPTracker* tracker = new SMPTracker(bline);

	// Perform tracking without errors
	tracker->Track(bunch1);

	Print("Final (w/o errors) ", bunch1);

	AcceleratorErrors err;
	ofstream logs("error.log");
	err.Report(&logs);
	// Transformations are applied in the order they are given
	// Order of rotations: 1. around x, 2. around y, 3. around z
	//
	// SetErrors: clear any existing local frame transformation
	// AddErrors: add to the existing offset/rotation
	//
	// ApplyShifts:    shifts in x,y,z as defined with SetErrors/AddErrors
	// ApplyRotations: rotations around x,y,z as defined with SetErrors/AddErrors
	//
	// pattern: type.name e.g.: Quadrupole.QFHE

	// random azimuthal error for all cavities
	//
	// arguments: rms_x,y,z, mean_x,y,z - default=0
	err.SetErrors(0, 100 * microradian);
	// arguments: beam line & pattern
	err.ApplyRotations(bline, "*.*CAV");

	// random transverse error for all cavities
	//
	err.AddErrors(300 * micrometer, 300 * micrometer);
	err.ApplyShifts(bline, "*.*CAV");

	// random transverse error for all quad
	// 200mu in blineFirst, 300mu in bline without blineFirst
	//
	err.SetErrors(300 * micrometer, 300 * micrometer);
	err.ApplyShifts(bline, "Quadrupole.*");
	//
	err.SetErrors(200 * micrometer, 200 * micrometer);
	err.ApplyShifts(blineFirst, "Quadrupole.*");
	//
	err.AddErrors(0, 0, 300 * microradian);
	err.ApplyRotations(blineFirst, "Quadrupole.*");

	// Perform tracking with errors
	tracker->Track(bunch2);

	Print("Final (with errors) ", bunch2);

	return 0;

}

// construct model and beam data from xttf-file fname
pair<AcceleratorModel*, BeamData*> ConstructModel(const string& fname)
{

	// bunch charge
	double qt = 2.0e+10;

	ofstream logs("construction.log");
	XTFFInterface mc(fname, qt, &logs);
	mc.ConstructGirders(true);

	pair<AcceleratorModel*, BeamData*> mb = mc.Parse();
	BeamData*         beam0 = mb.second;

	// The following quantities are not
	// int the TAPE file
	double gamma = beam0->p0 / MeV / ElectronMassMeV;
	beam0->emit_x = 8.0e-06 / gamma;
	beam0->emit_y = 0.02e-06 / gamma;
	beam0->charge = qt;
	beam0->sig_dp = 0.028;
	beam0->sig_z  = 300.0e-06;

	return mb;
}

// calculate emittances and print
void Print(string tag, SMPBunch* bunch)
{

	PSmoments S;
	bunch->GetMoments(S);
	double E = bunch->GetReferenceMomentum();
	E *= 1 + S.mean(ps_DP);
	double gam = E / MeV / ElectronMassMeV;
	double emitx = sqrt(S(0, 0) * S(1, 1) - pow(S(0, 1), 2));
	double emity = sqrt(S(2, 2) * S(3, 3) - pow(S(2, 3), 2));
	cout << tag << " energy " << E << endl;
	cout << tag << " gamma*emittance_x " << emitx * gam << endl;
	cout << tag << " gamma*emittance_y " << emity * gam << endl;

}
