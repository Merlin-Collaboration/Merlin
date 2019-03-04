/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "FrameModifier.h"

FrameModifier::FrameModifier(LatticeFrame* frame, const std::string& label) :
	LatticeFrame(label), subFrame(frame)
{
	LatticeFrame* oldSuperFrame = frame->SetSuperFrame(this);
	if(oldSuperFrame)
	{
		oldSuperFrame->ReplaceSubFrame(frame, this);
	}

	SetGeometry(frame->GetGeometry());
}

FrameModifier::~FrameModifier()
{
	Remove();
}

ModelElement* FrameModifier::Copy() const
{
	return subFrame->Copy();
}

const string& FrameModifier::GetType() const
{
	return subFrame->GetType();
}

void FrameModifier::Invalidate() const
{
	subFrame->Invalidate();
}

bool FrameModifier::IsBoundaryPlane(BoundaryPlane p, const LatticeFrame* aSubFrame) const
{
	assert(aSubFrame == subFrame);
	return true;
}

void FrameModifier::ConsolidateConstruction()
{
	subFrame->SetLocalPosition(0);
}

void FrameModifier::Remove()
{
	superFrame->ReplaceSubFrame(this, subFrame);
}
