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
// $Revision: 1.2.4.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_SMP_WakeFieldProcess
#define _h_SMP_WakeFieldProcess

#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamDynamics/SMPTracking/SMPBunchProcess.h"

#include <vector>
#include <typeinfo>

class WakePotentials;

namespace SMPTracking {

// Applies single-bunch long. and trans. wakefields to a SliceMacroParticles

class WakeFieldProcess : public SMPBunchProcess
{
public:

    enum ImpulseLocation {atCentre,atExit};

    WakeFieldProcess (int prio, double slice_width =1.0e-6,string aID="SMP WAKEFIELD");
    ~WakeFieldProcess();

    virtual void InitialiseProcess (Bunch& bunch);
    virtual void SetCurrentComponent (AcceleratorComponent& component);
    virtual void DoProcess (double ds);
    virtual double GetMaxAllowedStepSize () const;

    void ApplyImpulseAt(ImpulseLocation loc) { imploc=loc; }

    void IncludeTransverseWake(bool flg) { inc_tw = flg; }

private:

    void ApplyWakefield(double ds);

    ImpulseLocation imploc;
    double current_s;
    double impulse_s;
    double clen;
    const double dz; // slice width for binning

    bool inc_tw;
    void Init();
    void PrepLWake();
    void PrepSlices();

    std::vector<double> wake_z;
    std::vector<SMPBunch::iterator> sliceBoundaries;
    std::vector<double> slice_z;
    std::vector<double> slice_q;
    double bload;

    bool recalc;

    WakePotentials* currentWake;

};

}; // end namespace SMPTracking

#endif
