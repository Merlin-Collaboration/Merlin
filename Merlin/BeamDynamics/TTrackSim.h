/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/07 09:14:12 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_TTrackSim
#define _h_TTrackSim 1

#include "BeamDynamics/TrackingSimulation.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

// Transport process template function.
// Defines a transport BunchProcess which uses a concrete ComponentTracker
// class to track provide the tracking.
//
// Template argument T must be a valid concrete ComponentTracker class.

template<class T>
class TTrnsProc : public TBunchProc<__TYPENAME__ T::bunch_type> {
public:

    typedef __TYPENAME__ T::bunch_type bunch_type;
    typedef __TYPENAME__ T::B_Integrator integrator_type;
    typedef __TYPENAME__ T::ISetBase integrator_set_base;

    // Construct process using the default integrator set
    TTrnsProc()
            : TBunchProc<__TYPENAME__ T::bunch_type>("TRANSPORT",0), ctracker() {}

    // Construct process using an explicit integrator set
    TTrnsProc(const integrator_set_base& iset)
            : TBunchProc<__TYPENAME__ T::bunch_type>("TRANSPORT",0), ctracker(iset) {}

    // Register an addition (or override) integrator.
    bool RegisterIntegrator(integrator_type* intg) {
        return ctracker.Register(intg);
    }

    void InitialiseProcess (Bunch& bunch) {
        TBunchProc<bunch_type>::InitialiseProcess(bunch);
        if(this->currentBunch!=0) {
            this->active=true;
            ctracker.SetBunch(*this->currentBunch);
        }
    }

    void SetCurrentComponent (AcceleratorComponent& component) {
        if(this->active) {
            component.PrepareTracker(ctracker);
            this->currentComponent = &component;
        }
        else
            this->currentComponent=0;
    }

    void DoProcess(double ds) {
        if(this->active)
            ctracker.TrackStep(ds);
    }

    double GetMaxAllowedStepSize() const {
        return this->active ? ctracker.GetRemainingLength() : 0;
    }

private:

    T ctracker;
};

// template class to instantiate a TrackingSimulation class which
// uses a specific concrete ComponentTracker to track the
// associated bunch representation through a beamline.
//
// Template parameter T must be a valid concrete ComponentTracker class

template<class T>
class TTrackSim : public TrackingSimulation {
public:

    typedef __TYPENAME__ T::bunch_type bunch_type;
    typedef __TYPENAME__ bunch_type::particle_type particle_type;
    typedef TTrnsProc<T> transport_process;
    typedef __TYPENAME__ transport_process::integrator_type integrator_type;

    //	Constructor taking the beamline to be tracked and a
    //	pointer to the initial ParticleBunch. If bunch0=0
    //	(default) the initial beam must be specified by either a
    //	call to SetInitialBeam(), or to SetInitialBeamCtor (from
    //	TrackingSimulation).
    explicit TTrackSim (const AcceleratorModel::Beamline& bline,
                        bunch_type* bunch0 = 0, bool del = false);

    //	Constructor used for single particle tracking. The
    //	constructor takes the beamline to be tracked, and the
    //	initial particle.
    TTrackSim (const AcceleratorModel::Beamline& bline,
               const particle_type& p, double Pref);

    //	Constructor taking the beamline to be tracked and a
    //	pointer to the initial ParticleBunch. If bunch0=0
    //	(default) the initial beam must be specified by either a
    //	call to SetInitialBeam(), or to SetInitialBeamCtor (from
    //	TrackingSimulation).
    explicit TTrackSim (const AcceleratorModel::RingIterator& ring,
                        bunch_type* bunch0 = 0, bool del = false);

    //	Constructor used for single particle tracking. The
    //	constructor takes the beamline to be tracked, and the
    //	initial particle.
    TTrackSim (const AcceleratorModel::RingIterator& ring,
               const particle_type& p, double Pref);


    // Default constructor
    TTrackSim ();

    // Register an additional (or override) integrator.
    bool RegisterIntegrator(integrator_type* intg) {
        return transportProc->RegisterIntegrator(intg);
    }

    //	Sets the initial ParticleBunch for future tracking
    //	operations.
    void SetInitialBunch (bunch_type* pbunch0, bool del = false) {
        if(pbunch0!=0)
            SetInitialBunchCtor(MakeBunchCtor(pbunch0,del));
    }

    // Overrides the current bunch constructor and
    // tracks the supplied bunch.
    bunch_type* Track(bunch_type*);

    //	Sets the initial particle for single-particle tracking.
    void SetInitialParticle (const particle_type& p, double Pref);

    //	Returns a reference to the current tracked bunch.
    const bunch_type& GetTrackedBunch () const {
        return static_cast<const bunch_type&>(*bunch);
    }
    bunch_type& GetTrackedBunch () {
        return static_cast<bunch_type&>(*bunch);
    }

private:
    transport_process* transportProc;
};

// template implementation

template<class T>
TTrackSim<T>::TTrackSim (const AcceleratorModel::Beamline& bline,
                         bunch_type* bunch0, bool del)
        : TrackingSimulation(bline),transportProc(new transport_process())
{
    SetInitialBunch(bunch0,del);
    AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim (const AcceleratorModel::Beamline& bline,
                         const particle_type& p, double Pref)
        : TrackingSimulation(bline),transportProc(new transport_process())
{
    SetInitialParticle(p,Pref);
    AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim (const AcceleratorModel::RingIterator& ring,
                         bunch_type* bunch0, bool del)
        : TrackingSimulation(ring),transportProc(new transport_process())
{
    SetInitialBunch(bunch0,del);
    AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim (const AcceleratorModel::RingIterator& ring,
                         const particle_type& p, double Pref)
        : TrackingSimulation(ring),transportProc(new transport_process())
{
    SetInitialParticle(p,Pref);
    AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim ()
        : TrackingSimulation(),transportProc(new transport_process())
{
    AddProcess(transportProc);
}


template<class T>
void TTrackSim<T>::SetInitialParticle (const particle_type& p, double Pref)
{
    bunch_type* b = new bunch_type(Pref);
    b->AddParticle(p);
    SetInitialBunch(b,true);
}

template<class T>
__TYPENAME__ TTrackSim<T>::bunch_type* TTrackSim<T>::Track(bunch_type* aBunch)
{
    if(bunch!=0)
        delete bunch;

    bunch=aBunch;
    try {
        DoRun(false,true);
    }
    catch(...) {
        bunch=0;
        throw;
    }

    bunch=0;
    return aBunch;
}

#endif
