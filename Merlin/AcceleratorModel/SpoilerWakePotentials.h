/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
//  Created: June, 2006
// 
/////////////////////////////////////////////////////////////////////////

#ifndef SpoilerWakePotentials_h
#define SpoilerWakePotentials_h 1

#include "AcceleratorModel/WakePotentials.h"
#include "merlin_config.h"

//	Abstract class for calculating the longitudinal and
//	transverse single-bunch wakefield potentials (Greens
//	functions) with modes

class SpoilerWakePotentials : public WakePotentials
{
  

public:

SpoilerWakePotentials(int m) 
   : WakePotentials() { nmodes = m; }

    virtual ~SpoilerWakePotentials () {};
   
    virtual double Wlong (double s, int m) const = 0;
    virtual double Wtrans(double s, int m) const = 0;

protected:
   
    int nmodes;   
};

#endif
