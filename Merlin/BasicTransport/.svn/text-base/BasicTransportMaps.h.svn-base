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
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _H_BasicTransportMaps
#define _H_BasicTransportMaps 1

#include "BasicTransport/RTMap.h"

// BasicTransportMaps.h
//
// Contains declarations of global functions for constructing
// first- and second-order transport matrices. Both order maps
// are constructed and returned in a RTMap object.

RTMap* DriftTM(double s);
RTMap* SectorBendTM(double s, double h);
RTMap* SextupoleTM(double s, double K2);
RTMap* QuadrupoleTM(double s, double K1);
RTMap* GenSectorBendTM(double s, double h, double K1, double K2);
RTMap* PoleFaceTM(double h, double K1, double beta, double c, double fint, double hgap);

#endif
