/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ResponseMatrixGenerator.h"

ResponseMatrixGenerator::ResponseMatrixGenerator(Accelerator* acc1, const ROChannelArray& b, RWChannelArray& c, double
	eps1) :
	acc(acc1), bpms(b), cors(c), eps(eps1), data0(b.Size()), M(b.Size(), c.Size())
{
}

const RealMatrix& ResponseMatrixGenerator::GetMatrix() const
{
	return M;
}

const RealVector& ResponseMatrixGenerator::GetReference() const
{
	return data0;
}

const RealMatrix& ResponseMatrixGenerator::Generate(size_t ns)
{
	// generate reference
	acc->TrackBeam(ns);
	data0 = bpms;

	for(size_t icor = 0; icor < cors.Size(); icor++)
	{
		double defaultValue = cors.Read(icor);
		cors.Write(icor, defaultValue + eps);
		acc->TrackBeam(ns);
		cors.Write(icor, defaultValue);
		RealVector data(bpms);
		data -= data0;
		data /= eps;
		M.column(icor) = data;
	}

	return M;
}
