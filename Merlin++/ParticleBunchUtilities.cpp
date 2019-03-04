/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "merlin_config.h"
#include "ParticleBunchUtilities.h"

#include <vector>
#include <cmath>

#ifdef ENABLE_MPI
#include <mpi.h>
#include <fstream>
#include <iostream>
#endif

namespace
{

using namespace std;
using namespace ParticleTracking;

size_t TruncateZ(ParticleBunch& particles, double zmin, double zmax)
{
	size_t lost = 0;
	ParticleBunch::iterator p = particles.begin();
	while(p != particles.end())
	{
		double z = p->ct();
		if(z < zmin || z >= zmax)
		{
			p = particles.erase(p);
			lost++;
		}
		else
		{
			p++;
		}
	}
	return lost;
}

}

namespace ParticleTracking
{

// Sort the bunch in ascending z (ct) order, and return
// a vector of iterators which point to the equal-spaced
// bin boundaries defines by zmin to zmax in steps of dz
//
// Returns the number of particles removed from tails
// i.e. z<zmin || z>=zmax
//
// hdp contains the derivative of the distribution
// calculated using the Savitzky-Golay filter c
// If c is empty, then the derivative will be zero.
size_t ParticleBinList(ParticleBunch& bunch, double zmin, double zmax, size_t nbins,
	vector<ParticleBunch::iterator>& pbins, vector<double>& hd, vector<double>& hdp, vector<double>* c)
{
//cout << "In ParticleBinList" << endl;
//cout << zmin << "\t" << zmax << "\t" << nbins << endl;
	double dz = (zmax - zmin) / double(nbins);
	vector<ParticleBunch::iterator> bins;
	vector<double> hbins(nbins, 0);
	bins.reserve(nbins + 1);

	bunch.SortByCT();

	size_t lost = TruncateZ(bunch, zmin, zmax);

	ParticleBunch::iterator p = bunch.begin();
	bins.push_back(p);

	double z = zmin;
	double total = 0;
	size_t n;

	for(n = 0; n < nbins; n++)
	{
		z += dz;
		while(p != bunch.end() && p->ct() < z)
		{
			total++;
			hbins[n]++;
			p++;
		}
		bins.push_back(p);
	}

	if(p != bunch.end())
	{
#ifdef ENABLE_MPI
		cerr << "bad slicing in rank: " << MPI::COMM_WORLD.Get_rank() << endl;
#endif

#ifndef ENABLE_MPI
		cerr << "bad slicing" << endl;
#endif

		cerr << "z = " << z << " ct = " << p->ct() << " zmax = " << zmax << endl;
		//Dump out the bad bunch
		/*
		    ofstream* badbunch = new ofstream("badbunch.bunch");
		    bunch.Output(*badbunch);
		    badbunch->close();
		    delete badbunch;
		    cerr << "Output of the current bunch is to badbunch.bunch" << endl;
		 */
#ifndef ENABLE_MPI
		abort();
#endif

#ifdef ENABLE_MPI
		MPI::COMM_WORLD.Abort(1);
#endif
	}

	//	bins.push_back(p); // should be end()

	// normalise distribution
	// and apply filter

	vector<double> fbins(nbins, 0);
	vector<double> fpbins(nbins, 0);

	double a = 1 / total / dz;
	int w = c ? (c->size() - 1) / 2 : 0;
	size_t m;

	for(n = 0; n < nbins; n++)
	{
		fbins[n] = hbins[n] * a;
		if(c)
			//for(m=_MAX(0,int(n)-w); m<=_MIN(nbins,int(n)+w); m++)// ERROR! m can be set to nbins -> out of range!
			for(m = _MAX(0, int(n) - w); m < _MIN(nbins, size_t(n) + w); m++)   // This needs to be checked!
			{
				fpbins[n] += hbins[m] * (*c)[m - n + w] * a;
			}
	}

	pbins.swap(bins);
	hd.swap(fbins);
	hdp.swap(fpbins);

	return lost;
}

// Return the distribution of particles for the coordinate u.
// The distribution is returned as a binned histogram, with
// bin boundaries defined by umin to umax in steps of du.
// If truncate is true, then particles outside the range (umin,umax)
// are removed from bunch. If normalise is true, the histogram
// values are scaled to give a probability distribution with the
// property Sum{h_i*du} = 1.
//
size_t ParticleBunchDistribution(ParticleBunch& bunch, PScoord u, double umin, double umax, double du,
	vector<double>& bins, bool normalise, bool truncate)
{
	size_t nb = (umax - umin) / du;
	bins = vector<double>(nb, 0.0);
	size_t np0 = bunch.size();

	ParticleBunch::iterator p = bunch.begin();
	size_t lost = 0;
	double total = 0;

	while(p != bunch.end())
	{
		double uval = (*p)[u];
		if((uval < umin) || (uval > umax))
		{
			lost++;
			if(truncate)
			{
				p = bunch.erase(p);
				p--;
			}
		}
		else
		{
			size_t n = (uval - umin) / du;
			assert(n < nb);
			bins[n]++;
			total++;
		}
		p++;
	}
	if(normalise)
	{
		double factor = 1 / total / du;
		for(size_t i = 0; i < bins.size(); i++)
		{
			bins[i] *= factor;
		}
	}
	return lost;
}

}
