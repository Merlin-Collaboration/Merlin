/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include "MADInterface.h"
#include "LatticeFunctions.h"
#include "AcceleratorModel.h"
#include "BetatronTunes.h"
#include "ParticleBunchTypes.h"

#include "TransferMatrix.h"
#include "ClosedOrbit.h"

/* Read a TFS lattice with the MAD interface.
 *
 * Compare the tune found by the TransferMatrix and by FFT using BetatronTunes.
 *
 */

int main(int argc, char* argv[])
{

	const double beam_energy = 7000.0;

	// Find lattice file
	MADInterface* myMADinterface;
	string paths[] = {"../data/twiss.7.0tev.b1_new.tfs", "data/twiss.7.0tev.b1_new.tfs", "MerlinTests/data/twiss.7.0tev.b1_new.tfs"};

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

	// read lattice using MAD interface
	myMADinterface = new MADInterface(lattice_path, beam_energy);
	myMADinterface->TreatTypeAsDrift("RFCAVITY");
	AcceleratorModel* model = myMADinterface->ConstructModel();

	// find closed orbit
	double cscale = 1e-16;
	double delta = 1.0e-8;
	Particle p2(0);
	ClosedOrbit co(model, beam_energy);
	co.SetDelta(delta);
	co.ScaleBendPathLength(cscale);
	co.FindClosedOrbit(p2);

	// find transfer matrix
	RealMatrix M(6);
	TransferMatrix tm(model, beam_energy);
	tm.SetDelta(delta);
	tm.ScaleBendPathLength(cscale);
	tm.FindTM(M, p2);

	// calculate tune
	const double tm_qx = acos((M(0, 0) + M(1, 1)) / 2) / 2 / M_PI;
	const double tm_qy = acos((M(2, 2) + M(3, 3)) / 2) / 2 / M_PI;

	// measure tunes by fft
	BetatronTunes* tune =  new BetatronTunes(model, beam_energy);
	Particle p(0);
	p.x() = 1e-10;
	p.y() = 1e-10;

	tune->FindTunes(p, 512, false);
	const double fft_qx = tune->Qx, fft_qy = tune->Qy;

	// compare
	cout << "From TM" << endl;
	cout << tm_qx << " " << tm_qy << endl;
	cout << "From FFT" << endl;
	cout << fft_qx << " " << fft_qy << endl;

	double err_x = fabs(tm_qx - fft_qx), err_y = fabs(tm_qy - fft_qy);
	cout << "diffs " << err_x << " " << err_y << endl;

	assert(err_x < 1e-3);
	assert(err_y < 1e-2);

	delete myMADinterface;
	delete model;
	delete tune;

}
