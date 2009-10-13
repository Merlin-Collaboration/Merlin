/////////////////////////////////////////////////////////////////////////
// OneToOneCorrection
// Applies simply one-to-one steering of the BPMs using all the available
// corrections in the specified plane. The routine uses SVD and divides up
// the accelerator beamline in consecutive sections of nseg BPMs. The correction
// is applied 100% and is not iterated.
//
// Note that the current BeamDynamicsModel in the Accelerator is used for both
// the calcualtion of the response matrices and the estimate of the correction.
//
// The routine will attempt to use incremental tracking if possible.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_OneToOneCorrection
#define _h_OneToOneCorrection 1

#include "Accelerator.h"

void OneToOneCorrection(Accelerator* acc, Accelerator::Plane pxy, size_t nseg=20);

#endif

