/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ComponentFrame.h"

ComponentFrame::~ComponentFrame()
{
}

void ComponentFrame::Invalidate() const
{
}

const string& ComponentFrame::GetType() const
{
	_TYPESTR(ComponentFrame);
}

ModelElement* ComponentFrame::Copy() const
{
	return new ComponentFrame(*this);
}

bool ComponentFrame::IsBoundaryPlane(BoundaryPlane p, const LatticeFrame* aSubFrame) const
{
	// Should never be called!
	assert(false);
	return false;
}
