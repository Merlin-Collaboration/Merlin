/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ResponseMatrixGenerator
#define _h_ResponseMatrixGenerator

#include "TLAS.h"
#include "Channels.h"
#include "Accelerator.h"

// Utility class for constructing the response matrix for a given set of
// correctors and BPMs. Once constructed, the response matrix and reference
// trajectory can be calculated for any valid beam state.
class ResponseMatrixGenerator
{
public:

	ResponseMatrixGenerator(Accelerator* acc, const ROChannelArray& bpms, RWChannelArray& cors, double eps = 1.0e-06);

	// Generate the response matrix for the ns-th beam state
	const RealMatrix& Generate(size_t ns);

	const RealMatrix& GetMatrix() const;
	const RealVector& GetReference() const;

private:

	Accelerator* acc;
	const ROChannelArray& bpms;
	RWChannelArray& cors;
	double eps;

	RealVector data0;
	RealMatrix M;
};

#endif
