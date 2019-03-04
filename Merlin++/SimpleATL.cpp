/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <algorithm>
#include <iomanip>
#include "RandomNG.h"
#include "SimpleATL.h"

#undef ATL_XY

namespace
{

using namespace std;

inline void ResetSupport(AcceleratorSupport* s)
{
	s->Reset();
}

// function used to sort supports in acending arc position
inline bool AsZ(const AcceleratorSupport* s1, const AcceleratorSupport* s2)
{
	return s1->GetArcPosition() < s2->GetArcPosition();
}

struct ApplyATL
{

	ApplyATL(double A, double dt, vector<double>& gmy, double vibv, std::mt19937_64* rg) :
		AT(A * dt), y(0), z(0), yy(gmy.begin()), vv(vibv), rng(rg)
	{
	}

	void operator()(AcceleratorSupport* s)
	{
		// Perform the random walk
		double v = AT * (s->GetArcPosition() - z);
		y += normal_distribution<>{0, sqrt(v)} (*rng);
		*yy += y;

		// add random 'noise'
		double yv = vv != 0 ? normal_distribution<>{0, sqrt(vv)} (*rng) : 0.0;
		s->SetOffset(0, *yy + yv, 0);
		z = s->GetArcPosition();

		yy++;
	}

	// data members
	double AT;
	double y;
	double z;
	vector<double>::iterator yy;
	double vv;
	std::mt19937_64* rng;

};

struct DumpOffset
{

	DumpOffset(ostream& anOS) :
		os(anOS)
	{
	}

	void operator()(const AcceleratorSupport* s)
	{
		using std::setw;
		Vector3D offset = s->GetOffset();
		os << setw(12) << offset.x;
		os << setw(12) << offset.y;
		os << setw(12) << offset.z;
		os << endl;
	}

	ostream& os;

};

} // end of annonymous namespace

SimpleATL::SimpleATL(double anA, const AcceleratorSupportList& supports, double vrms) :
	t(0), A(anA), seed(0), vv(vrms * vrms), atlgm(supports.size(), 0.0), theSupports(supports)
{
	sort(theSupports.begin(), theSupports.end(), AsZ);
	rg = &RandomNG::getLocalGenerator(hash_string("ATL2D") + seed);
}

SimpleATL::~SimpleATL()
{
}

void SimpleATL::Reset()
{
	for_each(theSupports.begin(), theSupports.end(), ResetSupport);
	fill(atlgm.begin(), atlgm.end(), 0.0);
	t = 0;
}

double SimpleATL::DoStep(double dt)
{
	for_each(theSupports.begin(), theSupports.end(), ApplyATL(A, dt, atlgm, vv, rg));
	return t += dt;
}

void SimpleATL::RecordOffsets(std::ostream& os) const
{
	ios_base::fmtflags oldFlags = os.flags();
	int prec = os.precision();

	os.setf(ios_base::scientific, ios_base::floatfield);
	os.precision(4);

	for_each(theSupports.begin(), theSupports.end(), DumpOffset(os));

	os.flags(oldFlags);
	os.precision(prec);
}

double SimpleATL::GetTime() const
{
	return t;
}

void SimpleATL::SetRandomSeed(unsigned int nseed)
{
	seed = nseed;
	RandomNG::resetLocalGenerator(hash_string("ATL2D") + seed);
	rg = &RandomNG::getLocalGenerator(hash_string("ATL2D") + seed);
}

unsigned int SimpleATL::GetRandomSeed() const
{
	return seed;
}

void SimpleATL::ResetRandomSeed()
{
	RandomNG::resetLocalGenerator(hash_string("ATL2D") + seed);
	rg = &RandomNG::getLocalGenerator(hash_string("ATL2D") + seed);
}
