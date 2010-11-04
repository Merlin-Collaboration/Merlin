/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Created: June, 2006
/////////////////////////////////////////////////////////////////////////

#ifndef _h_SpoilerWakeProcess
#define _h_SpoilerWakeProcess

#include "merlin_config.h"
#include "Collimators/SpoilerWakePotentials.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "utility/StringPattern.h"
#include <vector>


//	Class for calculating the longitudinal and
//	transverse single-bunch wakefields 
//	for Spoilers with modes

namespace ParticleTracking {


class SpoilerWakeProcess : public WakeFieldProcess
{
public:

  SpoilerWakeProcess(int, int, size_t, double);
     ~SpoilerWakeProcess ();

    virtual void ApplyWakefield(double);
    virtual void CalculateWakeT (double, int);
    virtual void CalculateWakeL (double, int);

 private:

    double CalculateSm (int, int);
    double CalculateCm (int, int);
    
    int nmodes;
    
//    double Cm[5][1000];  
    double** Cm;
//    double Sm[5][1000];  
    double** Sm;
    double wake_s;
    double wake_c;
    
    double** wake_sl;
    double** wake_cl;
    double** wake_ct;
    double** wake_st;

//    ParticleBunch* currentBunch;
    SpoilerWakePotentials* spoiler_wake;  
};

}

#endif
