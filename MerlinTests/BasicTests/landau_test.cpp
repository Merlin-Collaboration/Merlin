/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include "RandomNG.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>

// This executable is called by landau_test.py
//
// The time (seconds since epoch) is used as the seed unless an int is
// passed as the first argument.


using namespace std;

int main(int argc, char* argv[])
{

	int seed = (int)time(nullptr);
	if (argc >=2)
	{
		seed = atoi(argv[1]);
	}

	cout << "Seed: " << seed << endl;

	RandomNG::init(seed);
	const int nthrows = 1e7;

	// if changing these, also change landau_test_gendata.py
	const size_t nbins = 2600;
	const double bin_min = -5, bin_max = 20;

	int hist[nbins+2] = {0}; // hist for Merlin
	// make some random numbers, and make the histogram
	for (int i = 0; i<nthrows; i++)
	{
		double x = RandomNG::landau();
		assert(!std::isnan(x));

		// beware, this can rollover when x is big
		int bin = ((x - bin_min) / (bin_max-bin_min) * (nbins)) +1; // +1 because bin zero for outliers
		// so handle end bins, by check against x, not bin
		if (x < bin_min)
		{
			bin = 0;
		}
		if (x > bin_max)
		{
			bin = nbins+1;
		}
		hist[bin] += 1;
	}

	// save histogram
	ofstream outfile("landau_test_output.dat");
	outfile << "#bin, count\n#From landau_test.cpp\n" ;
	for (size_t i = 0; i<nbins+2; i++)
	{
		outfile << i << " " << hist[i] <<endl;
	}

	return 0;
}


