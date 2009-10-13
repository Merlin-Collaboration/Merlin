/////////////////////////////////////////////////////////////////////////
// File ModelConstruction.h
// Routine to construct ILC Main Linac model
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_ModelConstruction
#define _h_ModelConstruction

#include <utility>
#include <string>

class AcceleratorModel;
class BeamData;

std::pair<AcceleratorModel*, BeamData*> ConstructModel(const std::string& fname, bool curved);

#endif

