/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_OneToOneCorrection
#define _h_OneToOneCorrection 1

#include "Accelerator.h"

// Applies simply one-to-one steering of the BPMs using all the available
// corrections in the specified plane. The routine uses SVD and divides up
// the accelerator beamline in consecutive sections of nseg BPMs. The correction
// is applied 100% and is not iterated.
//
// Note that the current BeamDynamicsModel in the Accelerator is used for both
// the calcualtion of the response matrices and the estimate of the correction.
//
// The routine will attempt to use incremental tracking if possible.
void OneToOneCorrection(Accelerator* acc, Accelerator::Plane pxy, size_t nseg = 20);

#endif
