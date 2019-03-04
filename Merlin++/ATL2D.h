/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ATL2D_h
#define ATL2D_h 1

#include "merlin_config.h"
#include <iostream>
#include "AcceleratorSupport.h"
#include "LinearAlgebra.h"

/**
 *
 *	Represents a 2D ATL model of ground motion. On each
 *	time step dT, ATL2D applies a random vertical (and
 *	horizontal if included) displacement to a set of
 *	AcceleratorSupports. The variance of the motion is
 *
 *	                 v = A.dT.L,
 *
 *	where A is the ATL constant, dT the time step and L the
 *	distance between supports.
 *
 *  ATL2D calculates correctly the correlation between the support point
 *  on a 2D plane, such that A.dT.L  holds for any two points, where L is
 *  the direct distance between those two points.
 *
 */

class ATL2D
{
public:

	enum ATLMode
	{
		increment,
		absolute

	};

	/**
	 *	Constructor taking the A constant and the list of
	 *	support structures.
	 */
	ATL2D(double anA, const AcceleratorSupportList& supports, const Point2D refPoint = Point2D(0, 0),
		ifstream* evecTFile = nullptr, ifstream* evalFile = nullptr);

	~ATL2D();

	/**
	 *	Reset the ground motion to zero Note this resets the
	 *	offset of all the AcceleratorSupports, and resets the
	 *	internal clock to zero.
	 */
	void Reset();

	/**
	 *	Perform a single step of dt seconds. Returns the current
	 *	simulated time.
	 *	@return Current simulated time (s)
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

	bool SetATLMode(const ATLMode mode);

	bool SetVibration(const double vrms);

	void RecordEigenSystem(ofstream* evecTFile, ofstream* evalFile);

private:

	double t;
	double A;
	unsigned seed;

	// Uncorrelated white-noise vibration variance
	double vv;

	AcceleratorSupportList theSupports;
	std::mt19937_64* rg;

	RealMatrix evecsT;
	RealVector evals;
	ATLMode atlMode;

	double Distance(const int n1, const int n2);
	double Distance(const int n1, const Point2D x2);
};

#endif
