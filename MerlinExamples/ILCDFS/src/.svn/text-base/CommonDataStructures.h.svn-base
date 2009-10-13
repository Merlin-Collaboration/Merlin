/////////////////////////////////////////////////////////////////////////
// CommonDataStructure.h
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

#ifndef _h_CommonDataStructures
#define _h_CommonDataStructures 1

#include "AcceleratorModel/ControlElements/Klystron.h"
#include "BeamModel/ReferenceParticle.h"
#include <iostream>
#include <vector>
#include <utility>

typedef std::vector<Klystron*> KlystronArray;
typedef std::vector<ReferenceParticle*> ReferenceParticleArray;
typedef std::pair<size_t,size_t> DFS_Segment;
typedef std::vector<size_t> IntegerArray;

inline std::ostream& operator<<(std::ostream& os, const DFS_Segment& seg)
{
	return os<<'['<<seg.first<<','<<seg.second<<']';
}

#endif
