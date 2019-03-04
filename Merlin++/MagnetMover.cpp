/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "MagnetMover.h"

const string& MagnetMover::GetType() const
{
	_TYPESTR(MagnetMover);
}

ModelElement* MagnetMover::Copy() const
{
	// Not sure what to do here!
	assert(false);
	return nullptr;
}

Transform3D MagnetMover::GetLocalFrameTransform() const
{
	return Transform3D(t) * LatticeFrame::GetLocalFrameTransform();
}
