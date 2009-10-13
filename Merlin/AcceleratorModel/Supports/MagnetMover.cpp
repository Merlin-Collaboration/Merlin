/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/Supports/MagnetMover.h"

const string& MagnetMover::GetType () const
{
    _TYPESTR(MagnetMover);
}

ModelElement* MagnetMover::Copy () const
{
    // Not sure what to do here!
    assert(false);
    return 0;
}

Transform3D MagnetMover::GetLocalFrameTransform () const
{
    return Transform3D(t)*LatticeFrame::GetLocalFrameTransform();
}

