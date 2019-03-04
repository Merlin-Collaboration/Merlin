/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _BunchConverter_h
#define _BunchConverter_h

#include "ParticleBunch.h"
#include "SMPBunch.h"

/*
 * Conversion functions between SMPBunch and ParticleBunch
 */

/**
 * Given an SPMBunch we construct a ParticleBunch with (on average) N particles.
 * The discrete ct, dp points in the original SMPBunch are smeared out in the ParticleBunch
 * with gauss(delta/2) where delta is given by the distance in either ct or dp.
 * adjust==true: force means(ParticleBunch)=means(SMPBunch)
 *
 */
ParticleTracking::ParticleBunch* SMPBunchConverter(SMPTracking::SMPBunch* SB, int N, bool adjust = true);

/**
 * Given a ParticleBunch we construct a SMPBunch with n_ct x n_dp particles
 * adjust==true: force means(ParticleBunch)=means(SMPBunch)
 * nSigZ,nSigDP number of sigmas to scan(integrate) the ParticleBunch, the sigmas are calculated from PB
 * nSigZ==0;nSigDP==0 scan from smallest to highest value in PB
 *
 */
SMPTracking::SMPBunch* ParticleBunchConverter(ParticleTracking::ParticleBunch* PB, int n_ct, int n_dp, double nSigZ =
	0, double nSigDP = 0, bool adjust = true);

#endif
