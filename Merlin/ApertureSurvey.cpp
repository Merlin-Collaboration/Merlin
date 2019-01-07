/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <fstream>

#include "ApertureSurvey.h"
#include "Collimator.h"
#include "ParticleBunch.h"

using namespace std;

static ofstream* open_file(const string file_name)
{
	ofstream* output_file = new ofstream(file_name.c_str());
	if(!output_file->good())
	{
		std::cerr << "Could not open ApertureSurvey file:" << file_name << std::endl;
		exit(EXIT_FAILURE);
	}
	return output_file;
}

void ApertureSurvey::ApertureSurvey(AcceleratorModel* model, string file_name, SurveyType mode, double step)
{
	// Simple wrapper around overloaded function that expects a ostream
	ofstream* output_file = open_file(file_name);
	ApertureSurvey(model, output_file, mode, step);
	delete output_file;
}

static void CheckAperture(Aperture* ap, double s, double *aps)
{
	// check the aperture a specific point in an element
	// uses bisection search

	const double delta = 1e-6; // search resolution
	const double max = 1.0;
	const double min = 0.0;

	// iterate through directions
	for(int dir = 0; dir < 4; dir++)
	{
		double xdir = 0, ydir = 0;
		if(dir == 0)
		{
			xdir = +1;
		}
		else if(dir == 1)
		{
			xdir = -1;
		}
		else if(dir == 2)
		{
			ydir = +1;
		}
		else if(dir == 3)
		{
			ydir = -1;
		}

		aps[dir] = 0;

		// scan for limit
		double below = min, above = max;

		while(above - below > delta)
		{
			double guess = (above + below) / 2;

			if(ap->CheckWithinApertureBoundaries(xdir * guess, ydir * guess, s))
			{
				below = guess;
			}
			else
			{
				above = guess;
			}
		}
		aps[dir] = (above + below) / 2;
	}
}

void ApertureSurvey::ApertureSurvey(AcceleratorModel* model, std::ostream* os, SurveyType mode, double step)
{
	(*os) << "#name\ttype\ts_end\tlength\ts\tap_px\tap_mx\tap_py\tap_my" << endl;

	// Some checks on the parameters
	size_t points = 0;
	double step_size = step;
	switch(mode)
	{
	case points_per_element:
		points = size_t(step);
		if(points < 1)
		{
			std::cerr << "ERROR: With mode=points_per_element, step must be at least 1" << std::endl;
			exit(EXIT_FAILURE);
		}
		break;
	case distance:
	case abs_distance:
		if(step_size <= 0)
		{
			std::cerr << "ERROR: step must be positive" << std::endl;
			exit(EXIT_FAILURE);
		}
		break;
	}

	double s = 0; // current element start
	double last_sample = 0 - step_size;
	double lims[4];

	// iterate though model
	for(AcceleratorModel::BeamlineIterator bi = model->GetBeamline().begin(); bi != model->GetBeamline().end(); bi++)
	{

		AcceleratorComponent *ac = &(*bi)->GetComponent();
		if(ac->GetLength() == 0)
		{
			continue;    // skip zero length elements
		}

		if(fabs(s - ac->GetComponentLatticePosition()) > 1e-6)
		{
			std::cerr << "ERROR: discrepancy between GetComponentLatticePosition() and sum of lengths" << std::endl;
			exit(EXIT_FAILURE);
		}

		// get aperture
		Aperture* ap =  ac->GetAperture();

		// find points within element to check, based on mode
		vector<double> zs;
		switch(mode)
		{
		case points_per_element:
			zs.push_back(0); // avoid divide by zero in loop if points=1
			for(size_t i = 1; i < points; i++)
			{
				zs.push_back(i * ac->GetLength() * (1.0 / (points - 1)));
			}

			break;
		case distance:
			for(double sample = 0; sample <= ac->GetLength(); sample += step_size)
			{
				zs.push_back(sample);
			}
			zs.push_back(ac->GetLength());
			break;
		case abs_distance:
			while(last_sample + step_size < s + ac->GetLength())
			{
				last_sample += step_size;
				zs.push_back(last_sample - s);
			}
			break;

		}

		// check aperture at each point in the element and output
		for(size_t zi = 0; zi < zs.size(); zi++)
		{
			double z = zs[zi];
			//cout << "call check_aperture(" << z << ")" << endl;
			if(ap != NULL)
			{
				CheckAperture(ap, z, lims);
			}
			else
			{
				lims[0] = lims[1] = lims[2] = lims[3] = 1;
			}

			(*os) << ac->GetName() << "\t";
			(*os) << ac->GetType() << "\t";
			(*os) << ac->GetComponentLatticePosition() + ac->GetLength() << "\t";
			(*os) << ac->GetLength() << "\t";
			(*os) << s + z << "\t";
			(*os) << lims[0] << "\t";
			(*os) << lims[1] << "\t";
			(*os) << lims[2] << "\t";
			(*os) << lims[3] << endl;
		}
		s += ac->GetLength();
	}
}
