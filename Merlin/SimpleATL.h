/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SimpleATL_h
#define SimpleATL_h 1

#include "merlin_config.h"
#include <iostream>
#include "AcceleratorSupport.h"

/**
 *	Represents a simple ATL model of ground motion. On each
 *	time step dT, SimpleATL applies a random vertical (and
 *	horizontal if included) displacement to a set of
 *	AcceleratorSupports. The variance of the motion is
 */
//	                 v = A.dT.L,
/**
 *  \f[
 *                   v = A \mathrm{d}T L
 *  \f]
 *	where A is the ATL constant, dT the time step and L the
 *	distance between supports.
 *
 *	Note that at present only the arc distance between
 *	successive supports  is used. For curved geometries this
 *	may introduce an error.
 */

class SimpleATL
{
public:
	/**
	 *	Constructor taking the A constant and the list of
	 *	support structures.
	 */
	SimpleATL(double anA, const AcceleratorSupportList& supports, double vrms = 0);

	~SimpleATL();

	/**
	 *	Reset the ground motion to zero Note this resets the
	 *	offset of all the AcceleratorSupports, and resets the
	 *	internal clock to zero.
	 */
	void Reset();

	/**
	 *	Perform a single step of dt seconds. Returns the current
	 *	simulated time.
	 */
	double DoStep(double dt);

	/**
	 *	Record the (x,y,z) offset of all the supports to the
	 *	specified stream.
	 */
	void RecordOffsets(std::ostream& os) const;

	/**
	 *	Returns the current simulated time (in seconds).
	 *	@return Current simulated time (s)
	 */
	double GetTime() const;

	/**
	 *	Sets the random seed to nseed.
	 */
	void SetRandomSeed(unsigned int nseed);

	/**
	 *	Returns the current random seed
	 *	@return Current random seed
	 */
	unsigned int GetRandomSeed() const;

	/**
	 *	Resets the random generator with the current random seed.
	 */
	void ResetRandomSeed();

private:

	double t;
	double A;
	unsigned seed;

	/**
	 * Uncorrelated white-noise vibration variance
	 */
	double vv;

	/**
	 * vector to store ATL ground motion
	 */
	std::vector<double> atlgm;

	AcceleratorSupportList theSupports;

	std::mt19937_64* rg;

	//Copy protection
	SimpleATL(const SimpleATL& rhs);
	SimpleATL& operator=(const SimpleATL& rhs);

};

#endif
