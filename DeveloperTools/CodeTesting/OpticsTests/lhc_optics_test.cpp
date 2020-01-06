/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "MADInterface.h"
#include "LatticeFunctions.h"
#include "DataTable.h"
#include "DataTableTFS.h"

using namespace std;

/* Read a TFS lattice with the MAD interface, measure the twiss parameters with
 * LatticeFunctionTable. Also read the data from the TFS file
 *
 * Check that the lattice has be read correctly, (names in the AcceleratorModel
 * match the names in the TFS file). Check that the values from LatticeFunctionTable
 * agree with those in the TFS file.
 *
 * Note that MADInterface deliberately ignores some zero length elements
 */

const vector<string> skip_types{"DRIFT", "RCOLLIMATOR", "INSTRUMENT", "TKICKER", "PLACEHOLDER", "VMONITOR", "HMONITOR",
								"KICKER", "MATRIX"};

bool check_diff(double x1, double x2, double da, double dr)
{
	if(x1 == 0 && x2 == 0)
	{
		return 1;
	}

	double absdiff = fabs(x1 - x2);
	double reldiff = fabs(x1 - x2) / max(x1, x2);
	if((absdiff > da) && (reldiff > dr))
	{
		cout << "absdiff " << absdiff << endl;
		cout << "reldiff " << reldiff << endl;
		return 0;
	}
	return 1;
}

int main(int argc, char* argv[])
{

	const double beam_energy = 7000.0;

	// Find lattice file
	MADInterface* myMADinterface;

	string lattice_path = find_data_file("twiss.7.0tev.b1_new.tfs");
	cout << "Lattice " << lattice_path << endl;

	// read lattice using MAD interface
	myMADinterface = new MADInterface(lattice_path, beam_energy);
	myMADinterface->TreatTypeAsDrift("RFCAVITY");
	AcceleratorModel* model = myMADinterface->ConstructModel();

	model->ReportModelStatistics(cout);

	//Calculate the main parameters (TWISS + dispersion) for ideal machine
	LatticeFunctionTable* twiss_table = new LatticeFunctionTable(model, beam_energy);

	twiss_table->AddFunction(1, 6, 3);
	twiss_table->AddFunction(2, 6, 3);
	twiss_table->AddFunction(3, 6, 3);
	twiss_table->AddFunction(4, 6, 3);
	twiss_table->AddFunction(6, 6, 3);

	auto mad_optics = make_unique<DataTable>(DataTableReaderTFS(lattice_path).Read());

	twiss_table->SetForceLongitudinalStability(true);
	twiss_table->Calculate();

	cout << setprecision(9);

	size_t x = 0; // counter for lattice table
	size_t xm = 0; // counter for mad tfs file
	for(const auto &bi : model->GetBeamline())
	{
		const AcceleratorComponent *ac = &(bi->GetComponent());
		// twiss_table has params at start of element
		// tfs has params at end, hence use so compare to previous

		if(x >= 1)
		{
			// Need to skip zero length DRIFTs and RCOLLIMATORs, as these get dropped by madinterface
			// so xm get extra increments
			while(true)
			{
				if(mad_optics->Get_d("L", xm) != 0)
					break;
				string mad_el_type = mad_optics->Get_s("KEYWORD", xm);
				if(mad_el_type == "MULTIPOLE" && mad_optics->Get_d("K0L", xm) == 0
					&& mad_optics->Get_d("K1L", xm) == 0
					&& mad_optics->Get_d("K2L", xm) == 0
					&& mad_optics->Get_d("K3L", xm) == 0
					&& mad_optics->Get_d("K4L", xm) == 0
					)
				{
					mad_el_type = "DRIFT";
				}

				if(none_of(skip_types.begin(), skip_types.end(), [&mad_el_type](const string st){
					return st == mad_el_type;
				}))
					break;
				xm++;
			}

			// check that we are on the correct element. (BPMs get there name manged, so ignore them)
			if(ac->GetName() != mad_optics->Get_s("NAME", xm) && mad_optics->Get_s("NAME", xm).substr(0, 3) != "BPM")
			{
				cout << "Names do not match at:" << x << " " << ac->GetName() << " != " << mad_optics->Get_s("NAME",
					x) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(1, 1, 1, x), mad_optics->Get_d("BETX", xm - 1), 5e-3, 1e-6))
			{
				cout << "Beta_x does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(1, 1, 1, x) << " != " << mad_optics->Get_d("BETX", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(3, 3, 2, x), mad_optics->Get_d("BETY", xm - 1), 5e-3, 1e-6))
			{
				cout << "Beta_y does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(3, 3, 2, x) << " != " << mad_optics->Get_d("BETY", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(-twiss_table->Value(1, 2, 1, x), mad_optics->Get_d("ALFX", xm - 1), 1e-3, 1e-6))
			{
				cout << "Alpha_x does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << -twiss_table->Value(1, 2, 1, x) << " != " << mad_optics->Get_d("ALFX", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(-twiss_table->Value(3, 4, 2, x), mad_optics->Get_d("ALFY", xm - 1), 1e-3, 1e-6))
			{
				cout << "Alpha_y does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << -twiss_table->Value(3, 4, 2, x) << " != " << mad_optics->Get_d("ALFY", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(1, 6, 3, x) / twiss_table->Value(6, 6, 3, x), mad_optics->Get_d("DX", xm
				- 1), 5e-3, 5e-2))
			{
				cout << "Disp_x does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(1, 6, 3, x) / twiss_table->Value(6, 6, 3, x) << " != "
					 << mad_optics->Get_d("DX", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(2, 6, 3, x) / twiss_table->Value(6, 6, 3, x), mad_optics->Get_d("DPX",
				xm - 1), 2e-3, 1e-2))
			{
				cout << "Disp_px does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(2, 6, 3, x) / twiss_table->Value(6, 6, 3, x) << " != "
					 << mad_optics->Get_d("DPX", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(3, 6, 3, x) / twiss_table->Value(6, 6, 3, x), mad_optics->Get_d("DY", xm
				- 1), 5e-3, 5e-2))
			{
				cout << "Disp_y does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(3, 6, 3, x) / twiss_table->Value(6, 6, 3, x) << " != "
					 << mad_optics->Get_d("DY", xm - 1) << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(4, 6, 3, x) / twiss_table->Value(6, 6, 3, x), mad_optics->Get_d("DPY",
				xm - 1), 2e-3, 1e-2))
			{
				cout << "Disp_py does not match at:" << x << " " << mad_optics->Get_s("NAME", xm)
					 << " " << twiss_table->Value(4, 6, 3, x) / twiss_table->Value(6, 6, 3, x) << " != "
					 << mad_optics->Get_d("DPY", xm - 1) << endl;
				exit(1);
			}
		}
		x++;
		xm++;
	}

	delete myMADinterface;
	delete model;
	delete twiss_table;
}
