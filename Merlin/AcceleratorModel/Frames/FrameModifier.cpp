/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

// FrameModifier
#include "AcceleratorModel/Frames/FrameModifier.h"

FrameModifier::FrameModifier (LatticeFrame* frame, const std::string& label)
        : LatticeFrame(label),subFrame(frame)
{
    LatticeFrame* oldSuperFrame = frame->SetSuperFrame(this);
    if(oldSuperFrame)
        oldSuperFrame->ReplaceSubFrame(frame,this);

    SetGeometry(frame->GetGeometry());
}

FrameModifier::~FrameModifier()
{
    Remove();
}

ModelElement* FrameModifier::Copy () const
{
    return subFrame->Copy();
}

const string& FrameModifier::GetType () const
{
    return subFrame->GetType();
}

void FrameModifier::Invalidate () const
{
    subFrame->Invalidate();
}

bool FrameModifier::IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const
{
    assert(aSubFrame==subFrame);
    return true;
}

void FrameModifier::ConsolidateConstruction ()
{
    subFrame->SetLocalPosition(0);
}

void FrameModifier::Remove()
{
    superFrame->ReplaceSubFrame(this,subFrame);
}
