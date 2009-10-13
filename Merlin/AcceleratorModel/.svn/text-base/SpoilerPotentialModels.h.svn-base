/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3.2 (2008)
// 
// Copyright: see Merlin/copyright.txt
//
// Created: DK 25.2.2008
//    see BeamDynamics/ParticleTracking/SpoilerWakeProcess.cpp
/////////////////////////////////////////////////////////////////////////

#ifndef SpoilerPotentialModels_h
#define SpoilerPotentialModels_h 1

#include "merlin_config.h"
#include "AcceleratorModel/SpoilerWakePotentials.h"

//----------------------------------------------------------------------------------------------------------
//   The geometric wake potential:
//   Steeply tapered collimator moving from aperture b to aperture a 
//   Ref.:
//   B.W.Zotter and S.A.Kheifets, Impedances and Wakes in High-Energy Particle Accelerators, 
//   World Scientific (1998)
//----------------------------------------------------------------------------------------------------------
class TaperedCollimatorPotentials: public SpoilerWakePotentials  {
public:

    TaperedCollimatorPotentials(int m, double aa, double bb);
    ~TaperedCollimatorPotentials();
    virtual double Wlong (double z, int m) const;
    virtual double Wtrans (double z, int m) const;
private:
    double* coeff;
    double a, b; 
};
     
//------------------------------------------------------------------------------------------------------------
//      the resistive wake potentials  (in MKS ssytem)
//
//------------------------------------------------------------------------------------------------------------
class ResistiveWakePotentials: public SpoilerWakePotentials {
public:
    ResistiveWakePotentials(int m, double r, double s, double l) ;
    ~ResistiveWakePotentials();
    virtual double Wlong (double z, int m) const;
    virtual double Wtrans (double z, int m) const;
   
private:
   double* coeff;
   double rad, sigma, length;
};
#endif
