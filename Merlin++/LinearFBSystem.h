/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef LinearFBSystem_h
#define LinearFBSystem_h 1

#include "merlin_config.h"

#include "TLAS.h"
#include <queue>
#include <list>
#include "Channels.h"
#include "LinearAlgebra.h"

using namespace TLAS;

/**
 *	A simple linear feedback correction algorithm. On each
 *	application, the actuator channels (A) are incremented
 *	using the following linear equation:
 *
 *	A = A-g*(Mi*(S-S0))
 *
 *	where g is the gain, S are the current signal values and
 *	S0 are the desired signal values (set points). Mi is a
 *	pseudo-inverse matrix of the response matrix M defined by
 *
 *	S=M*A
 *
 *	Mi is calculated using SVD.
 */

class LinearFBSystem
{
public:

	LinearFBSystem(std::vector<ROChannel*>& sigs, std::vector<RWChannel*>& acts, const RealMatrix& M);
	LinearFBSystem(ROChannelArray& sigs, RWChannelArray& acts, const RealMatrix& M);

	~LinearFBSystem();

	void SignalsToSetpoints();
	void StoreActuators() const;
	void RestoreActuators();
	void SetResponseMatrix(const RealMatrix& M);
	void SetGain(double g);
	double GetGain() const;
	void Apply();
	void SetSetpoints(const RealVector& S0);
	double GetActuatorRMS() const;
	double GetSignalRMS() const;
	int GetNumSignals() const;
	int GetNumActuators() const;
	void SetPulseDelay(int n);

private:

	double gain;
	ROChannelArray signals;
	RWChannelArray actuators;
	RealVector setpoints;

	LinearFBSystem(const LinearFBSystem &right);
	const LinearFBSystem & operator=(const LinearFBSystem &right);

	mutable RealVector* cached_actuators;

	SVDMatrix<double>* Mi;
	// to allow for possible actuator pulse delays, we use a queue
	mutable std::queue<RealVector>* actuatorQueue;
};

#endif
