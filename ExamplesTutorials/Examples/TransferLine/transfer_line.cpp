/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <fstream>
#include <memory>

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

#include "ParticleDistributionGenerator.h"
#include "AcceleratorModelConstructor.h"
#include "StandardMultipoles.h"
#include "SectorBend.h"
#include "Drift.h"
#include "BeamData.h"

#include "ClosedOrbit.h"
#include "LatticeFunctions.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;

int main()
{
	RandomNG::init(5);

	// Create ring cell
	AcceleratorModelConstructor ring_constructor;

	double quad_k = 0.0081;
	double quad_l = 5;
	double bend_angle = 2 * pi / 200;
	double bend_h = bend_angle / quad_l;
	double bend_l = 5;

	double beam_energy = 20 * GeV;
	double beam_momentum = sqrt(pow(beam_energy, 2) - pow(ProtonMassMeV * MeV, 2));
	double rigid = beam_energy / eV / SpeedOfLight;

	auto qf = new Quadrupole("QF", quad_l, quad_k * rigid);
	auto qd = new Quadrupole("QD", quad_l, -quad_k * rigid);
	auto mb = new SectorBend("MB", bend_l, bend_h, rigid * bend_h);
	auto d = new Drift("D", 3);

	ring_constructor.AppendComponent(qf);
	ring_constructor.AppendComponent(d);
	ring_constructor.AppendComponent(mb);
	ring_constructor.AppendComponent(d);
	ring_constructor.AppendComponent(qd);
	ring_constructor.AppendComponent(qd);
	ring_constructor.AppendComponent(d);
	ring_constructor.AppendComponent(mb);
	ring_constructor.AppendComponent(d);
	ring_constructor.AppendComponent(qf);

	AcceleratorModel* ring_lattice = ring_constructor.GetModel();

	// Calculate beta and dispersion functions
	LatticeFunctionTable latticeFunctions = LatticeFunctionTable(ring_lattice, beam_momentum);
	latticeFunctions.SetForceLongitudinalStability(true);
	latticeFunctions.Calculate();
	ofstream latticeFunctionLog("example_transferline.out");
	latticeFunctions.PrintTable(latticeFunctionLog);

	double beta_x_ring = latticeFunctions.Value(1, 1, 1, 0);
	double beta_y_ring = latticeFunctions.Value(3, 3, 2, 0);

	cout << "Ring beta_x " << beta_x_ring << " beta_y " << beta_y_ring << endl;

	// Create transfer line
	AcceleratorModelConstructor transfer_constructor;

	double quadt_l = 4;
	vector<double> quadt_k = {1.26827e-02,  -1.44910e-02, -1.53154e-02, 1.58761e-02};
	//vector<double> quadt_k = {1.0e-02,  -1.0e-02, -1.0e-02, 1.0e-02};
	auto q1 = new Quadrupole("Q1", quadt_l, quadt_k[0] * rigid);
	auto q2 = new Quadrupole("Q2", quadt_l, quadt_k[1] * rigid);
	auto q3 = new Quadrupole("Q3", quadt_l, quadt_k[2] * rigid);
	auto q4 = new Quadrupole("Q4", quadt_l, quadt_k[3] * rigid);
	auto dt = new Drift("DT", 3);

	transfer_constructor.AppendComponent(dt);
	transfer_constructor.AppendComponent(q1);
	transfer_constructor.AppendComponent(dt);
	transfer_constructor.AppendComponent(q2);
	transfer_constructor.AppendComponent(dt);
	transfer_constructor.AppendComponent(q3);
	transfer_constructor.AppendComponent(dt);
	transfer_constructor.AppendComponent(q4);
	transfer_constructor.AppendComponent(dt);

	AcceleratorModel* transfer_lattice = transfer_constructor.GetModel();

	// create bunch
	double emittance = 1e-6;
	int npart = 10000;
	BeamData mybeam;
	mybeam.p0 = beam_momentum;
	mybeam.beta_x = 30;
	mybeam.beta_y = 20;
	mybeam.emit_x = emittance;
	mybeam.emit_y = emittance;

	ParticleBunch myBunch(npart, NormalParticleDistributionGenerator(), mybeam);

	ofstream bunchout0("example_transferline_bunch_init.out");
	myBunch.Output(bunchout0, true);

	//track
	auto bline_ring = ring_lattice->GetRing();
	auto bline_transfer = transfer_lattice->GetBeamline();

	ParticleTracker tracker_ring(bline_ring, &myBunch);
	ParticleTracker tracker_transfer(bline_transfer, &myBunch);

	cout << "Bunch size init: " << myBunch.size() << endl;
	tracker_transfer.Track(&myBunch);

	cout << "Bunch size inject: " << myBunch.size() << endl;
	ofstream bunchout1("example_transferline_bunch_inject.out");
	myBunch.Output(bunchout1, true);

	for(int turn = 0; turn < 100; turn++)
	{
		tracker_ring.Track(&myBunch);
	}
	cout << "Bunch size end: " << myBunch.size() << endl;
	ofstream bunchout2("example_transferline_bunch_end.out");
	myBunch.Output(bunchout2, true);

	delete ring_lattice;
	delete transfer_lattice;
	return 0;
}
