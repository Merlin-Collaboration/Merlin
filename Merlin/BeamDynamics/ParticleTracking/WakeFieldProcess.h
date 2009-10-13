/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2008/03/03 13:58:26 $
// $Revision: 1.6.4.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_WakeFieldProcess
#define _h_WakeFieldProcess

#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "AcceleratorModel/WakePotentials.h"
#include "utility/StringPattern.h"
#include <vector>
#include <typeinfo>

class WakePotentials;

namespace ParticleTracking {

// class WakefieldProcess
// Implements single-bunch wakefields.
//
// Modified by A.Wolski 12/2/2003
// to include the derivative of the charge distribution
// for handling CSR wake fields.
//
// Modified by D.Kruecker 18.2.2008
// to be used as base class for other wakefield types (spoiler,coupler,...)

class WakeFieldProcess : public ParticleBunchProcess {
public:

    enum ImpulseLocation {atCentre,atExit};

    WakeFieldProcess (int prio, size_t nb =100, double ns = 3.0, string aID = "WAKEFIELD");

    ~WakeFieldProcess();

    virtual void InitialiseProcess (Bunch& bunch);
    virtual void SetCurrentComponent (AcceleratorComponent& component);
    virtual void DoProcess (double ds);
    virtual double GetMaxAllowedStepSize () const;

    void ApplyImpulseAt(ImpulseLocation loc) { imploc=loc; }
    void IncludeTransverseWake(bool flg) { inc_tw = flg; }
    void DumpSliceCentroids(ostream&) const;
    void SetFilter(int n, int m, int d);

    //private: 
    protected:

    ImpulseLocation imploc;
    double current_s;
    double impulse_s;
    double clen;
    size_t nbins;
    double nsig;
    bool inc_tw;

    void Init();
    size_t CalculateQdist();
    virtual void CalculateWakeL();
    virtual void CalculateWakeT();
    virtual void ApplyWakefield(double ds);

    WakePotentials* currentWake;

    std::vector<ParticleBunch::iterator> bunchSlices;
    std::vector<double> Qd;
    std::vector<double> Qdp;
    std::vector<double>* filter;
    std::vector<double> wake_z;
    std::vector<double> wake_x;
    std::vector<double> wake_y;

    double zmin,zmax,dz;
    bool recalc;
	size_t oldBunchLen;
};

void savgol(vector<double>& c, int nl, int nr, int ld, int m);

}; // end namespace ParticleTracking

#endif
