/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2019 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Components.h"
#include "AcceleratorModelConstructor.h"
#include "SymplecticIntegrators.h"
#include "MADInterface.h"

#include "ParticleTracker.h"
#include "ParticleBunchTypes.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

const std::string whitespace = " \t\f\v\n\r";

std::vector<std::string> split_line(std::string line)
{
	// Read line into vector of strings
	// split on whitespace unless quoted
	std::vector<std::string> sl;
	std::string current;
	current.reserve(64);
	bool double_quote_on = false;

	for(auto c : line)
	{
		if(c == '"')
		{
			double_quote_on = !double_quote_on;
		}
		else if(!double_quote_on && whitespace.find(c) != std::string::npos)
		{
			if(current.length())
			{
				sl.push_back(current);
				current =  "";
			}
		}
		else
		{
			current += c;
		}
	}
	if(current.length())
	{
		sl.push_back(current);
	}

	return sl;
}

int main(int argc, char* argv[])
{
	//Beam energy (GeV) 7000,3500,450 etc
	//double beam_energy = 7000.0;

	if(argc < 2)
	{
		cerr << "Usage:" << argv[0] << " infile" << endl;
		return 1;
	}
	string infile_name(argv[1]);

	std::ifstream in_file(infile_name);
	std::string line;
	map<string, string> settings;
	while(getline(in_file, line))
	{
		line = line.substr(0, line.find("#"));
		vector<string> words = split_line(line);
		if(words.size() == 0)
		{
			continue;
		}
		if(words[0] == "set")
		{
			settings[words.at(1)] = words.at(2);
		}
	}
	double beam_energy = stod(settings["energy"]);
	const double brho = beam_energy / eV / SpeedOfLight;

	if(settings["particle"] != "proton")
	{
		cerr << "merlin_track.cpp requires particle = proton" << endl;
		return 1;
	}

	ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1);
	AcceleratorComponent* element = nullptr;

	bool use_madinterface = false;
	string madinterface_tfs_path;

	in_file.clear();
	in_file.seekg(0);
	while(getline(in_file, line))
	{
		//cout << line << endl;
		line = line.substr(0, line.find("#"));
		vector<string> words = split_line(line);
		if(words.size() == 0)
		{
			continue;
		}
		if(words[0] == "particle")
		{
			Particle p(0);
			p.x() = stod(words.at(1));
			p.xp() = stod(words.at(2));
			p.y() = stod(words.at(3));
			p.yp() = stod(words.at(4));
			p.ct() = stod(words.at(5));
			p.dp() = stod(words.at(6));
			myBunch->push_back(p);
		}
		else if(words[0] == "set")
		{
			continue; // handled above
		}
		else if(words[0] == "madinterface")
		{
			use_madinterface = true;
			madinterface_tfs_path = words.at(1);
		}
		else if(words[0] == "drift")
		{
			element = new Drift("d1", stod(words.at(1)));
		}
		else if(words[0] == "quad")
		{
			double len = stod(words.at(1));
			double k1 = stod(words.at(2)) * brho;
			element = new Quadrupole("q1", len, k1);
		}
		else if(words[0] == "sbend")
		{
			double len = stod(words.at(1));
			double angle = stod(words.at(2));
			double h = angle / len;
			double b0 = brho * h;
			double k1l = stod(words.at(3));
			double k2l = stod(words.at(4));
			auto element_sb = new SectorBend("b1", len, h, b0);
			cout << "k1=" << k1l << " k2=" << k2l << endl;
			if(k1l)
				element_sb->SetB1(brho * k1l);

			//if(k2l)
			//	element_sb->SetBn(2, brho * k2l);
			element = element_sb;
		}
		else if(words[0] == "vkick")
		{
			double len = stod(words.at(1));
			double kick = brho * stod(words.at(2));
			if(len != 0)
				kick /= len;
			element = new YCor("b1", len, kick);
		}
		else if(words[0] == "hkick")
		{
			double len = stod(words.at(1));
			double kick = brho * stod(words.at(2));
			if(len != 0)
				kick /= len;
			element = new XCor("b1", len, -kick);
		}
		else
		{
			cerr << "Unknown line: \"" << line << "\"" << endl;
			return 1;
		}
	}

	AcceleratorModel* model;
	if(use_madinterface)
	{
		if(element != nullptr)
		{
			cerr << "Found element and madinterface keyword in " << infile_name << endl;
			return 1;
		}

		MADInterface* myMADinterface = new MADInterface(madinterface_tfs_path, beam_energy);
		model = myMADinterface->ConstructModel();
		delete myMADinterface;
	}
	else
	{
		if(element == nullptr)
		{
			cerr << "No element in " << infile_name << endl;
			return 1;
		}

		AcceleratorModelConstructor* construct = new AcceleratorModelConstructor();
		construct->NewModel();

		construct->AppendComponent(element);

		model = construct->GetModel();
		delete construct;
	}

	/*********************************************************************
	 *	PARTICLE TRACKER
	 *********************************************************************/

	AcceleratorModel::RingIterator bline = model->GetRing();
	ParticleTracker* tracker = new ParticleTracker(bline, myBunch);

	if(settings["int_mode"] == "")
	{
		cout << "Using Integrator: SYMPLECTIC (default)" << endl;
		tracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
	}
	else if(settings["int_mode"] == "symplectic")
	{
		cout << "Using Integrator: SYMPLECTIC" << endl;
		tracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
	}
	else if(settings["int_mode"] == "symplectic_ef")
	{
		cout << "Using Integrator: SYMPLECTIC" << endl;
		tracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
		tracker->RegisterIntegrator(new SYMPLECTIC::SectorBendCI_ef);
	}
	else if(settings["int_mode"] == "transport")
	{
		cout << "Using Integrator: TRANSPORT" << endl;
		tracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
	}
	else
	{
		cout << "Unknown Integrator: " << settings["int_mode"] << endl;
		return 1;
	}

	cout << "Tracking" << endl;
	tracker->Track(myBunch);

	string outbunch_fname = "merlin_out.dat";
	ofstream outbunch(outbunch_fname);

	if(!outbunch.good())
	{
		cerr << "Could not open " << outbunch_fname << endl;
		exit(1);
	}
	outbunch << "#p px y py ct pt\n" << endl;
	outbunch.precision(10);
	for(auto &&p : *myBunch)
	{
		outbunch << p.x() << " "
				 << p.xp() << " "
				 << p.y() << " "
				 << p.yp() << " "
				 << p.ct() << " "
				 << p.dp()
				 << endl;
	}
	outbunch.close();
	delete tracker;
	delete model;

	return 0;
}
