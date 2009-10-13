/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/15 13:43:32 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////


#include "MADInterface/ConstructSrot.h"
#include "AcceleratorModel/Frames/PatchFrame.h"

ComponentFrame* ConstructSrot(double angle, const std::string& name)
{
    GeometryPatch* gp = new GeometryPatch;
    gp->RotateZ(angle);
    return new PatchFrame(gp,name);
}

ComponentFrame* ConstructXrot(double angle, const std::string& name)
{
    if(angle) {
        GeometryPatch* gp = new GeometryPatch;
        gp->RotateX(angle);
        return new PatchFrame(gp,name);
    }
    else
        return new PatchFrame(0,name);
}
