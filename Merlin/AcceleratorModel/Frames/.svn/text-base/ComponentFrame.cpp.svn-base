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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

// ComponentFrame
#include "AcceleratorModel/Frames/ComponentFrame.h"

ComponentFrame::~ComponentFrame ()
{}

void ComponentFrame::Invalidate () const
{}

const string& ComponentFrame::GetType () const
{
    _TYPESTR(ComponentFrame);
}

ModelElement* ComponentFrame::Copy () const
{
    return new ComponentFrame(*this);
}

bool ComponentFrame::IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const
{
    // Should never be called!
    assert(false);
    return false;
}
