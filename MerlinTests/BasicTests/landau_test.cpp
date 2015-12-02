#include "../tests.h"
#include "Random/RandomNG.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

// This test should fail occasionally (10-20% of runs).
// It histograms the result of calling the RandomNG::landau() function
// and compares the result to an expected pdf calculated using GSL by
// the script landau_test_gendata.py
// The number of bins deviating by more than N sigma are recorded
// and if the are significantly more that expected (by a normal
// distribution), the test fails.
//
// If the test fails run it a few times. If it consistently fails then
// there is an actual issue. Bins out by more than 3 sigma are
// displayed, repeated failure of the same bin would be suspicious.
//
// The time (seconds since epoch) is used as the seed unless an int is
// passed as the first argument.


using namespace std;

int main(int argc, char* argv[])
{

	int seed = (int)time(NULL);
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
	double l_dist[nbins+2] = {0}; // pdf from gsl

	// Test could be run from number of directories
	string paths[] = {"../data/landau_test_dist.dat", "data/landau_test_dist.dat", "MerlinTests/data/landau_test_dist.dat"};
	ifstream dist_file;
	for (size_t i=0; i<3; i++){
		dist_file.open(paths[i].c_str());
		if (dist_file){break;}
	}
	if (! dist_file)
	{	
		cout << "Failed to open dist file"<< endl;
		exit(1);
	}
	
	// Load expected values	
	size_t dist_bin;
	double dist_val;
	string line;
	while (getline(dist_file, line)){
		if (line == "") continue;
		if (line[0] == '#')continue;
		istringstream liness(line);
		liness >> dist_bin >> dist_val;
		//cout << "@" << dist_bin << " " << dist_val << endl;
		if (dist_bin >= (nbins+2)){
			cout << "Invalid bin ("<< dist_bin<<") in landau_test_dist.dat";
			exit(1);
		}
		l_dist[dist_bin] = dist_val*nthrows;
	}

	// make some random numbers, and make the histogram	
	for (int i = 0; i<nthrows; i++){
		double x = RandomNG::landau();
		assert(!isnan(x));
		
		// beware, this can rollover when x is big
		int bin = ((x - bin_min) / (bin_max-bin_min) * (nbins)) +1; // +1 because bin zero for outliers
		// so handle end bins, by check against x, not bin
		if (x < bin_min) bin = 0;
		if (x > bin_max) bin = nbins+1;
		hist[bin] += 1;

	}

	// Calculate sigma from sqrt(<n>)
	// Count number of bins out by more than n-sigma
	const int ns = 4;
	// https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule, could use erf() in c++11
	const double tails[] = {0, 0.682689492137086, 0.954499736103642, 0.997300203936740,
	                           0.999936657516334, 0.999999426696856, 0.999999998026825};
	int over_ns[ns+1] = {0};

	for (size_t i = 0; i<nbins+2; i++){
		// should use poisson stats or ks test, otherwise bins with small counts can give bad results
		// Currently assume gaussian, don't allow sigma to be below 1.0
		double exp_err = max(sqrt(l_dist[i]), 1.0);

		double diff = hist[i] - l_dist[i];
		//cout << i << " " << hist[i] << " " << l_dist[i] << endl;
		//cout << i << " " << hist[i] << " " << l_dist[i] << "+-" << exp_err << endl;
		//cout << i << " " <<  fabs(diff)/exp_err << endl;
		
		for(int j=1; j<ns+1; j++){
			if (fabs(diff)/exp_err > j) {
				over_ns[j]++;
				// display the worst bins
				if (j>=3) {
					cout << endl;
					cout << "big error # bin, merlin, expected, error/sigma" << endl;
					cout << i << " " << hist[i] << " " << l_dist[i] << "+-" << exp_err << " " << diff/exp_err << endl;
				}
			}
		}
	}

	// put limit on fraction of bins out by more than n-sigma
	// bad if 30% more than expected
	// allow some leeway for low count bins, by using over_ns[j]-1)/(nbins+2)
	
	bool bad = 0;
	double tol = 1.3;
	cout << endl << "Bins out by more than";
	for(int j=1; j<ns+1; j++){
		cout << endl << "    "<< j <<" sigma: "<< over_ns[j] << " ("<< (double)over_ns[j]/(nbins+2)*100  << " %)";
		if ((double)over_ns[j]/(nbins+2) > (1-tails[j])*tol) {cout << "** bad **"; bad=1;}
	}
	cout << endl;
	
	return bad;
}


