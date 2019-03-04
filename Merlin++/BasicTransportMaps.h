/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _H_BasicTransportMaps
#define _H_BasicTransportMaps 1

#include "RTMap.h"

/*
 * BasicTransportMaps.h
 *
 * Contains declarations of global functions for constructing
 * first- and second-order transport matrices. Both order maps
 * are constructed and returned in a RTMap object.
 */

RTMap* DriftTM(double s);
RTMap* SectorBendTM(double s, double h);
RTMap* SextupoleTM(double s, double K2);
RTMap* QuadrupoleTM(double s, double K1);
RTMap* GenSectorBendTM(double s, double h, double K1, double K2);
RTMap* PoleFaceTM(double h, double K1, double beta, double c, double fint, double hgap, bool ent);

#endif
