/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SMPBunchConstructor_h
#define SMPBunchConstructor_h 1

#include "merlin_config.h"
#include "BeamData.h"
#include "BunchConstructor.h"
#include "SMPBunch.h"

namespace SMPTracking
{

class SMPBunchConstructor: public BunchConstructor
{
public:

	/**
	 * Constructor taking the beam definition, the number of z slices
	 * (ns) and the number of SMPs per slice. Total number of SMPs
	 * generated will be nsm*ns
	 */
	SMPBunchConstructor(const BeamData& beam, size_t ns, size_t nsm);

	~SMPBunchConstructor();

	virtual Bunch* ConstructBunch(int bunchIndex = 0) const;

	/**
	 *	Returns typed particle bunch.
	 */
	SMPBunch* ConstructSMPBunch() const;

private:

	size_t ns, np;
	BeamData beamdat;
	double nSigZ, nSigDP;
};

} // end namespace SMPTracking

#endif
