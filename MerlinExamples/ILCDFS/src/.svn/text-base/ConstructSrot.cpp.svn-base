/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
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
