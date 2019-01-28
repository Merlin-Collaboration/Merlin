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

/* Read a TFS lattice with the MAD interface, measure the twiss parameters with
 * LatticeFunctionTable. Also read the data from the TFS file
 *
 * Check that the lattice has be read correctly, (names in the AcceleratorModel
 * match the names in the TFS file). Check that the values from LatticeFunctionTable
 * agree with those in the TFS file.
 *
 * Note that MADInterface deliberately ignores some zero length elements
 */

struct tfs_line
{
	string name;
	string keyword;
	double s, l;
	double betx, bety, alfx, alfy;
	double k0l, k1l, k2l, k3l, k4l;

};

inline string StripQuotes(const string& text)
{
	return text.substr(1, text.length() - 2);
}

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

	model->ReportModelStatistics(cout);

	//Calculate the main parameters (TWISS + dispersion) for ideal machine
	LatticeFunctionTable* twiss_table = new LatticeFunctionTable(model, beam_energy);

	twiss_table->AddFunction(1, 6, 3);
	twiss_table->AddFunction(2, 6, 3);
	twiss_table->AddFunction(3, 6, 3);
	twiss_table->AddFunction(4, 6, 3);
	twiss_table->AddFunction(6, 6, 3);

	// read twiss parameters from the lattice file
	vector<tfs_line> tfs_table;
	ifstream tfs_file(lattice_path.c_str());
	string line;
	string dummy;
	while(getline(tfs_file, line))
	{
		while(line.length() > 0 && line[0] == ' ')
		{
			line = line.substr(1, line.length() - 1);
		}
		if(line == "")
		{
			continue;
		}
		if(line[0] == '@' || line[0] == '*' || line[0] == '$')
		{
			continue;
		}
		istringstream liness(line);
		tfs_line this_line;
		//(0, 'NAME'), (1, 'KEYWORD'), (2, 'S'), (3, 'L'), (4, 'KS'), (5, 'KSL'), (6, 'K0L'), (7, 'K1L'), (8, 'K2L'), (9, 'K3L'), (10, 'K4L'), (11, 'K1S'), (12, 'K2S'), (13, 'K3S'), (14, 'K4S'), (15, 'HKICK'), (16, 'VKICK'), (17, 'BETX'), (18, 'BETY'), (19, 'ALFX'), (20, 'ALFY'), (21, 'MUX'), (22, 'MUY'), (23, 'DX'), (24, 'DY'), (25, 'DPX'), (26, 'DPY'), (27, 'R11'), (28, 'R12'), (29, 'R22'), (30, 'R21'), (31, 'X'), (32, 'PX'), (33, 'Y'), (34, 'PY'), (35, 'T'), (36, 'DELTAP'), (37, 'VOLT'), (38, 'LAG'), (39, 'HARMON'), (40, 'FREQ'), (41, 'E1'), (42, 'E2'), (43, 'APERTYPE'), (44, 'APER_1'), (45, 'APER_2'), (46, 'APER_3'), (47, 'APER_4'), (48, 'TILT'), (49, 'ANGLE'), (50, 'ASSEMBLY_ID'), (51, 'MECH_SEP')
		for(int i = 0; i < 52; i++)
		{
			switch(i)
			{
			case 0:
				liness >> this_line.name;
				this_line.name = StripQuotes(this_line.name);
				break;
			case 1:
				liness >> this_line.keyword;
				this_line.keyword = StripQuotes(this_line.keyword);
				break;
			case 2:
				liness >> this_line.s;
				break;
			case 3:
				liness >> this_line.l;
				break;
			case 6:
				liness >> this_line.k0l;
				break;
			case 7:
				liness >> this_line.k1l;
				break;
			case 8:
				liness >> this_line.k2l;
				break;
			case 9:
				liness >> this_line.k3l;
				break;
			case 10:
				liness >> this_line.k4l;
				break;
			case 17:
				liness >> this_line.betx;
				break;
			case 18:
				liness >> this_line.bety;
				break;
			case 19:
				liness >> this_line.alfx;
				break;
			case 20:
				liness >> this_line.alfy;
				break;
			default:
				liness >> dummy;
				break;

			}
		}
		//liness >> this_line.name >> this_line.keyword;
		// Mimic MAD interface ignoring certain elements
		if(this_line.keyword == "MULTIPOLE")
		{
			if(this_line.k0l != 0)
			{
				this_line.keyword = "SBEND";
			}
			else if(this_line.k1l != 0)
			{
				this_line.keyword = "QUADRUPOLE";
			}
			else if(this_line.k2l != 0)
			{
				this_line.keyword = "SEXTUPOLE";
			}
			else if(this_line.k3l != 0)
			{
				this_line.keyword = "OCTUPOLE";
			}
			else if(this_line.k4l != 0)
			{
				this_line.keyword = "DECAPOLE";
			}
			else
			{
				this_line.keyword = "DRIFT";
			}
		}
		if(this_line.keyword == "INSTRUMENT")
		{
			this_line.keyword = "DRIFT";
		}
		if(this_line.keyword == "TKICKER")
		{
			this_line.keyword = "DRIFT";
		}
		if(this_line.keyword == "PLACEHOLDER")
		{
			this_line.keyword = "DRIFT";
		}

		if((this_line.keyword == "RCOLLIMATOR" ||
			this_line.keyword == "DRIFT")
			&& this_line.l == 0)
		{
			continue;
		}

		tfs_table.push_back(this_line);
	}

	double bscale1 = 1e-22;

	//loop
	while(true)
	{
		cout << "Trying bscale: " << bscale1 << "\t";
		cout.flush();
		twiss_table->ScaleBendPathLength(bscale1);
		cout << "Calculating" << endl;
		twiss_table->Calculate();
		cout << "Done calculating" << endl;
		if(!std::isnan(twiss_table->Value(1, 1, 1, 0)))
		{
			cout << "Success!" << endl;
			break;
		}

		bscale1 *= 2;
		cout << "Fail :(" << endl;
		if(bscale1 > 1e-18)
		{
			cout << "Giving up" << endl;
			exit(1);
		}
	}

	/*AddFunction(1,0,0); // closed orbit: x
	 *     AddFunction(2,0,0); // closed orbit: px
	 *     AddFunction(3,0,0); // closed orbit: y
	 *     AddFunction(4,0,0); // closed orbit: py
	 *     AddFunction(5,0,0); // closed orbit: ct
	 *     AddFunction(6,0,0); // closed orbit: dp
	 *     AddFunction(1,1,1); // beta_x
	 *     AddFunction(1,2,1); // -alfa_x
	 *     AddFunction(3,3,2); // beta_y
	 *     AddFunction(3,4,2); // -alfa_y
	 *                                     */

	cout << setprecision(9);

	size_t x = 0;
	//double worst_betx = 0, worst_bety = 0, worst_alfx = 0, worst_alfy = 0;
	//int worst_betx_n = 0, worst_bety_n = 0, worst_alfx_n = 0, worst_alfy_n = 0;
	for(AcceleratorModel::BeamlineIterator bi = model->GetBeamline().begin(); bi != model->GetBeamline().end(); bi++)
	{
		AcceleratorComponent *ac = &(*bi)->GetComponent();
		// twiss_table has params at start of element
		// tfs has params at end, hence use tfs_table[x-1]
		if(x >= 1)
		{
			//cout << ac->GetName() << "\t" ;//<< twiss_table->Value(1,1,1,x) << " ";
			//cout << tfs_table[x].name <<endl;//" " << tfs_table[x-1].betx << endl;

			// check that we are on the correct element. (BPMs get there name manged, so ignore them)
			if(ac->GetName() != tfs_table[x].name && tfs_table[x].name.substr(0, 3) != "BPM")
			{
				cout << "Names do not match at:" << x << " " << ac->GetName() << " != " << tfs_table[x].name << endl;
				exit(1);
			}

			//cout << ac->GetName() << "\t" << twiss_table->Value(1,2,1,x) << " " << tfs_table[x-1].alfx << endl;
			// check params

			if(!check_diff(twiss_table->Value(1, 1, 1, x), tfs_table[x - 1].betx, 5e-3, 1e-6))
			{
				cout << "Beta_x does not match at:" << x << " " << tfs_table[x].name
					 << " " << twiss_table->Value(1, 1, 1, x) << " != " << tfs_table[x - 1].betx << endl;
				exit(1);
			}

			if(!check_diff(twiss_table->Value(3, 3, 2, x), tfs_table[x - 1].bety, 1e-3, 1e-6))
			{
				cout << "Beta_y does not match at:" << x << " " << tfs_table[x].name
					 << " " << twiss_table->Value(3, 3, 2, x) << " != " << tfs_table[x - 1].bety << endl;
				exit(1);
			}

			if(!check_diff(-twiss_table->Value(1, 2, 1, x), tfs_table[x - 1].alfx, 1e-3, 1e-6))
			{
				cout << "Alpha_x does not match at:" << x << " " << tfs_table[x].name
					 << " " << -twiss_table->Value(1, 2, 1, x) << " != " << tfs_table[x - 1].alfx << endl;
				exit(1);
			}

			if(!check_diff(-twiss_table->Value(3, 4, 2, x), tfs_table[x - 1].alfy, 1e-3, 1e-6))
			{
				cout << "Alpha_y does not match at:" << x << " " << tfs_table[x].name
					 << " " << -twiss_table->Value(3, 4, 2, x) << " != " << tfs_table[x - 1].alfy << endl;
				exit(1);
			}

		}
		x++;
	}

	delete myMADinterface;
	delete model;
	delete twiss_table;

}
