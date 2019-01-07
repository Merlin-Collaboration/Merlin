/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "merlin_config.h"
#include "AcceleratorModel.h"
#include "XTFFInterface.h"
#include "BeamData.h"
#include "PhysicalConstants.h"
#include "TWRFStructure.h"
#include "Aperture.h"
#include "TeslaWakePotential.h"
#include "PatchFrame.h"
#include <vector>

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

pair<AcceleratorModel*, BeamData*> ConstructModel(const string& fname)
{
	double qt = 2.0e+10; // single-bunch charge

	ofstream logs("construction.log");
	XTFFInterface mc(fname, qt, &logs);
	mc.ConstructGirders(true);

//	mc.TreatTypeAsDrift("LCAV");

	pair<AcceleratorModel*, BeamData*> mb = mc.Parse();
	BeamData* beam0 = mb.second;
	AcceleratorModel* model = mb.first;

	// The following quantities are not
	// int the TAPE file

	double gamma = beam0->p0 / MeV / ElectronMassMeV;
	beam0->emit_x = 8.0e-06 / gamma;
	beam0->emit_y = 0.02e-06 / gamma;
	beam0->charge = qt == 0 ? 1.0 : qt;
	beam0->sig_dp = 0.028;
	beam0->sig_z  = 300.0e-06;

	// Now add wakefields to cavities

	vector<TWRFStructure*> cavities;
	model->ExtractTypedElements(cavities, "CAV*"); // linac cavities only
	TeslaWakePotentials* wake = new TeslaWakePotentials;

	Aperture* iris = new CircularAperture(0.035); // TESLA iris aperture

	for(vector<TWRFStructure*>::iterator c = cavities.begin(); c != cavities.end(); c++)
	{
		(*c)->SetWakePotentials(wake);
		(*c)->SetAperture(iris);
	}

	return mb;
}
