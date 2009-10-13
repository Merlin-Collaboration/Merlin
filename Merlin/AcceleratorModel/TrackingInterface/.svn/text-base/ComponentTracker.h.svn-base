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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ComponentTracker_h
#define ComponentTracker_h 1

#include "merlin_config.h"
#include "AcceleratorModel/TrackingInterface/ComponentIntegrator.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include <map>

//  class ComponentTracker
//
//	ComponentTracker provides the primary interface for tracking
//	operations through a series of AcceleratorComponents. An object of
//	class ComponentTracker can be thought of as a collection of
//	Integrator objects, each of which is responsible for
//	"tracking" some beam-like representation through a
//	specific component. Selection of the correct Integrator
//	is performed via the callback mechanism SelectIntegrator(),
//	which is called directly by Accelerator
//	AcceleratorComponent::PrepareTracker.
//  Concrete Tracker objects
//	should supply the (concrete) Integrators, which when
//	identified, are then passed back down to the concrete
//	class via a call the virtual function SetCurrent
//	Integrator.
//
//	The callback mechanism used by Tracker replaces the
//	canonical virtual function mechanism. It is a form of
//	the visitor pattern, but without the circular dependency
//	which is normally associated with that pattern.  Thus
//	the accelerator model representation can be easily
//	extended without the need to add new visitor functions.

class ComponentTracker
{
public:

    //	Tracker state
    typedef enum {undefined, initialised, tracking, finished} State;

    //	Exception thrown when no valid integrator is identified.
    struct UnknownComponent {};

    //	Destructor
    virtual ~ComponentTracker ();

    //	Track the entire current accelerator component. Initial
    //	tracker state must be initialised, and final state is
    //	finished.
    void Track ();

    //	Track a step ds through the current component, and
    //	return the remaining distance to the exit boundary.
    //	Initial state must be either initialised or tracking.
    //	Final state is either tracking or finished. If ds tracks
    //	beyond the current AcceleratorComponent geometry, a
    //	BeyondExtent exception is thrown.
    double TrackStep (double ds);

    //	Returns the current state of the Tracker.
    ComponentTracker::State GetState () const;

    //	Reset the Tracker. The state is automatically set to
    //	undefined.
    void Reset ();

    //	Returns the remaining tracking distance. Current state
    //	must be initialised, tracking or finished.
    double GetRemainingLength () const;

    //	Returns the total integrated length. State must be
    //	initialised, tracking or finished.
    double GetIntegratedLength () const;

    //	Function called by AcceleratorComponent objects to
    //	select the correct Integrator for component. Returns
    //	true if an Integrator object is found for index,
    //	otherwise false.
    bool SelectIntegrator (int index, AcceleratorComponent& component);

    //	Function operator overload. Tracks the specified
    //	AcceleratorComponent in one step.
    void operator () (AcceleratorComponent* component);

    // IntegratorSet class
    class IntegratorSet {
    public:
        typedef std::map< int,ComponentIntegrator* > IMap;
        ~IntegratorSet ();

        // Adds the specified integrator. Returns true
        // if the the integrator replaces an existing one,
        // otherwise false.
        bool Add (ComponentIntegrator* intg);

        // Return integrator n (or NULL if there is
        // no associate integrator).
        ComponentIntegrator* GetIntegrator (int n);

    private:
        IMap itsMap;
    };

protected:

    //	Constructor(s) made protected to prevent instantiation
    ComponentTracker ();
    explicit ComponentTracker (IntegratorSet* iset);

    //	Protected virtual function called by SelectIntegrator()
    //	when it has successfully located a valid integrator.
    //	Concrete Tracker objects should override this function
    //	to perform any initialisation appropriate on the
    //	(concrete) integrator.
    virtual void InitialiseIntegrator (ComponentIntegrator*);

    //	Used to register an integrator. Returns true if this
    //	integrator overrides an existing one, otherwise false.
    //	Register() should be called by the constructors of
    //	concrete Tracker classes in order to setup the required
    //	integrators.
    bool Register (ComponentIntegrator* intg);

private:

    // Dissable copy construction
    ComponentTracker(const ComponentTracker&);

    State itsState;
    ComponentIntegrator* integrator; // current integrator

    IntegratorSet* iSet;
};

// template class TBunchCMPTracker
//
// This template class can be used to instantiate concrete
// ComponentTracker classes which are used to track a specific
// Bunch object.

template<class _B>
class TBunchCMPTracker : public ComponentTracker {
public:

    typedef _B bunch_type;

    // Bunch specific integrators
class B_Integrator : public ComponentIntegrator {
    public:
        void SetBunch(_B& aBunch) { currentBunch = &aBunch; }
        virtual double Track(double ds) {
            double s = ComponentIntegrator::Track(ds);
            currentBunch->IncrReferenceTime(ds);
            return s;
        }
    protected:
        _B* currentBunch;
    };

    // Template for component-specific bunch integrator
    // concrete BunchType integrators should inherit from
    // an instantiation of this template.
    template<class _C>
class Integrator : public B_Integrator {
    public:
        typedef _C ComponentType;
    protected:
        // virtual function override
        void SetCurrentComponent (AcceleratorComponent& c) {
            ComponentIntegrator::SetCurrentComponent(c);
            currentComponent=static_cast<_C*>(&c);
        }

        int GetComponentIndex () const {
            return ComponentType::ID;
        }

        _C* currentComponent;
    };

    // Methods

    // Registration of bunch integrators
    bool Register(B_Integrator* ci) {
        return ComponentTracker::Register(ci);
    }

    // Set the bunch to be tracked.
    void SetBunch(_B& aBunch) {
        currentBunch = &aBunch;
    }

    // Integrator set definition.
    class ISetBase {
    public:
        virtual void Init(TBunchCMPTracker&) const =0;
    };

    // Default constructor (uses default integrator set)
    TBunchCMPTracker() : ComponentTracker(), currentBunch(0) {
        defIS->Init(*this);
    }

    // Constructor taking explicit integrator set
    explicit TBunchCMPTracker(const ISetBase& iset)
            : ComponentTracker(), currentBunch(0) {
        iset.Init(*this);
    }

    static void SetDefaultIntegratorSet(ISetBase* iset) {
        if(defIS!=0)
            delete defIS;
        defIS=iset;
    }

protected:

    static ISetBase* defIS;


    // virtual function override
    void InitialiseIntegrator (ComponentIntegrator* ci) {
        ComponentTracker::InitialiseIntegrator(ci);
        static_cast<B_Integrator*>(ci)->SetBunch(*currentBunch);
    }

    _B* currentBunch;
};

// macros for constructing integrator sets
#define DECL_INTG_SET(T,S) struct S : public T::ISetBase { void Init(T& ct) const; }; 
#define DEF_INTG_SET(T,S) void S::Init(T& ct) const {
#define ADD_INTG(iname) ct.Register(new iname ());
#define END_INTG_SET };
#define MAKE_DEF_INTG_SET(T,S) T::ISetBase* T::defIS = new S ();

// implementations

inline ComponentTracker::State ComponentTracker::GetState () const
{
    return itsState;
}

inline void ComponentTracker::operator () (AcceleratorComponent* component)
{
    component->PrepareTracker(*this);
    Track();
}

// Utility macros for tracker operations
#define _PREPTRACK(trk,S) if(!trk.SelectIntegrator(ID,*this)) S::PrepareTracker(trk); 


#endif
