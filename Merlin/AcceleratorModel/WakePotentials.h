/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2008/05/22 21:28:17 $
// $Revision: 1.3.4.3 $
// 
/////////////////////////////////////////////////////////////////////////
//
// Modified by D.Kruecker 18.2.2008
// to be used as base class for other wake potentials
// see SpoilerWakeProcess 

#ifndef WakePotentials_h
#define WakePotentials_h 1

#include "merlin_config.h"
#include <string>
#include <iostream>

#include "BeamDynamics/BunchProcess.h"

using namespace std;

//	Abstract class for calculating the longitudinal and
//	transverse single-bunch wakefield potentials (Greens
//	functions).

class WakePotentials
{
public:

    WakePotentials() : expectedProcess(0), csr(false){}
    virtual ~WakePotentials () {};

//    virtual double Wlong (double z)  const = 0;
//    virtual double Wtrans (double z) const = 0;
    virtual double Wlong (double z) const { 
    	cout<<"WakePotentials::Wlong: Virtual function called!"<<endl; 
    	abort();
    	return 0;
    };
    virtual double Wtrans (double z) const { 
    	cout<<"WakePotentials::Wtrans: Virtual function called!"<<endl; 
    	return 0;
    };
    bool Is_CSR () const {
        return csr;
    }
    BunchProcess* GetExpectedProcess(){return expectedProcess;};
    void          SetExpectedProcess(BunchProcess* p){expectedProcess=p;};
        
protected:
    bool csr;
    BunchProcess* expectedProcess;
};

#endif
