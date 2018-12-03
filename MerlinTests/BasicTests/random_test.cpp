/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "../tests.h"
#include "RandomNG.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <chrono>
#include <unordered_map>
#include "DataTable.h"
#include "DataTableTFS.h"

// This executable is called by random_test.py
//
// The time (seconds since epoch) is used as the seed unless an int is
// passed as the first argument.

using namespace std;

enum dist_names
{
	normal,
	normalc,
	uniform,
	poisson,
	landau

};
const unordered_map<string, dist_names> dist_name_map {
	{"normal", normal},
	{"normalc", normalc},
	{"uniform", uniform},
	{"poisson", poisson},
	{"landau", landau},
};

int main(int argc, char* argv[])
{

	int seed = 0;
	if(argc >= 2)
	{
		seed = atoi(argv[1]);
	}

	if(seed == 0)
	{
		seed = (int) time(nullptr);
	}

	RandomNG::init(seed);

	const int nthrows =  1e6;

	cout << "Seed: ";
	for(auto s : RandomNG::getSeed())
	{
		cout << s << " ";
	}
	cout << endl;

	std::unique_ptr<DataTable> random_setup(DataTableReaderTFS(find_data_file("random_test_setup.tfs")).Read());

	DataTable results;

	for(const auto row : *random_setup)
	{
		results.AddColumn(row.Get_s("name"), 'i');
	}

	// if changing these, also change landau_test_gendata.py
	const size_t nbins = 2000;

	for(size_t i = 0; i < nbins + 2; i++)
	{
		results.AddRow();
	}

	for(const auto row : *random_setup)
	{
		const string name = row.Get_s("name");
		const string dist = row.Get_s("dist");
		const auto dist_type = dist_name_map.at(dist);
		const double var1 = row.Get_d("var1");
		const double var2 = row.Get_d("var2");
		const double var3 = row.Get_d("var3");
		const double bin_min = row.Get_d("range_min");
		const double bin_max = row.Get_d("range_max");
		double x;
		bool is_cont = true;

		vector<double> start_vals;

		auto start = chrono::steady_clock::now();
		for(int i = 0; i < nthrows; i++)
		{
			switch(dist_type)
			{
			case normal:
				x = RandomNG::normal(var1, var2);
				break;
			case normalc:
				x = RandomNG::normal(var1, var2, var3);
				break;
			case uniform:
				x = RandomNG::uniform(var1, var2);
				break;
			case poisson:
				is_cont = false;
				x = RandomNG::poisson(var1);
				break;
			case landau:
				x = RandomNG::landau();
				break;
			default:
				cerr << "Unknown distribution: " << dist << endl;
				exit(1);
			}

			assert(!std::isnan(x));

			if(i < 10)
			{
				start_vals.push_back(x);
			}

			int bin;
			if(is_cont)
			{
				// so handle end bins, by check against x, not bin
				if(x < bin_min)
				{
					bin = 0;
				}
				else if(x > bin_max)
				{
					bin = nbins + 1;
				}
				else
				{
					bin = ((x - bin_min) / (bin_max - bin_min) * (nbins)) + 1; // +1 because bin zero for outliers
				}
			}
			else
			{
				// For a discrete distribution we use integer bins
				// and truncate at the bin max
				int bin_max_d = bin_max;
				if(x < bin_min)
				{
					bin = 0;
				}
				else if(x >= bin_max_d)
				{
					bin = bin_max_d - bin_min + 1;
				}
				else
				{
					bin = x - bin_min + 1;
				}
			}
			results.Set(name, bin, results.Get_i(name, bin) + 1);
		}
		auto end = chrono::steady_clock::now();
		chrono::duration<double, nano> diff = end - start;
		cout << "Time per value from  " << name << ": " << diff.count() / nthrows << " ns" << endl;
		for(auto x : start_vals)
		{
			cout << x << " ";
		}
		cout << endl;
	}

	ofstream outfile("random_test_output.dat");
	outfile << "# bins=" << nbins << " throws=" << nthrows << " seed=" << seed << endl;
	results.OutputAscii(outfile);

	return 0;
}
