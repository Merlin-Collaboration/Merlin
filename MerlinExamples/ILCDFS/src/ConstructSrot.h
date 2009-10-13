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

#ifndef _h_ConstructSrot
#define _h_ConstructSrot 1

#include <string>
#include "AcceleratorModel/Frames/ComponentFrame.h"

ComponentFrame* ConstructSrot(double angle, const std::string& name);
ComponentFrame* ConstructXrot(double angle, const std::string& name);

#endif