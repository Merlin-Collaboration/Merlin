/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "SMPComponentTracker.h"
#include "SMPStdIntegrators.h"

namespace SMPTracking
{

DEF_INTG_SET(SMPComponentTracker, StdISet)
ADD_INTG(DriftCI)
ADD_INTG(SectorBendCI)
ADD_INTG(RectMultipoleCI)
ADD_INTG(TWRFStructureCI)
ADD_INTG(MonitorCI)
ADD_INTG(SolenoidCI)
ADD_INTG(MarkerCI)
END_INTG_SET

} //end namespace SMPTracking

MAKE_DEF_INTG_SET(SMPTracking::SMPComponentTracker, SMPTracking::StdISet)
