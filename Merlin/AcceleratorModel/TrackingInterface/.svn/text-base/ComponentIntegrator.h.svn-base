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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ComponentIntegrator_h
#define ComponentIntegrator_h 1

#include "merlin_config.h"
#include "NumericalUtils/utils.h"
#include "NumericalUtils/Range.h"
#include "AcceleratorModel/AcceleratorComponent.h"

//  Class ComponentIntegrator
//
//	Abstract class used by ComponentTracker to track through an
//	AcceleratorComponent. Integrators encapsulate the
//	various algorithms used for tracking.

class ComponentIntegrator
{
public:

    ComponentIntegrator ();
    virtual ~ComponentIntegrator ();

    //	Track the entire component.
    void TrackAll ();

    //	Tracks a single step ds through the current component.
    //  Returns the remaining length.
    virtual double Track (double ds);

    //	Returns the component index for this integrator.
    virtual int GetComponentIndex () const = 0;

    //	Returns the total integrated length through the current
    //	component.
    double GetIntegratedLength () const;

    //	Returns the remaining length to integrate.
    double GetRemainingLength () const;

    //	Returns true if the step ds is valid (allowed)
    bool IsValidStep (double ds) const;

    //  Returns the remaining valid step range for this integrator.
    //  Default behaviour returns [0,GetRemainingLength()]
    virtual FloatRange GetValidStepRange() const;

    //	Component boundary checks
    //  Returns true if GetIntegratedLength()==0
    bool AtEntrance() const;

    //  Returns true if step == GetRemainingLength()
    bool AtExit(double step=0) const;

protected:

    friend class ComponentTracker;

    //  Sets the current component to track
    //  Can only be called by ComponentTracker via friendship
    virtual void SetCurrentComponent (AcceleratorComponent& aComponent);

    //  Perform the tracking for step ds. Concrete integrators
    //  must supply this function.
    virtual void TrackStep(double ds) =0;

    // functions for applying entrance and exit field maps.
    virtual void TrackEntrance() {};
    virtual void TrackExit() {};

private:

    double cur_S; // current integrated length
    double tot_S; // total length

    //  The current component
    AcceleratorComponent* component;
};

inline ComponentIntegrator::ComponentIntegrator ()
        : cur_S(0),component(0)
{}

inline ComponentIntegrator::~ComponentIntegrator ()
{}

inline void ComponentIntegrator::TrackAll ()
{
    Track(component->GetLength());
}

inline void ComponentIntegrator::SetCurrentComponent (AcceleratorComponent& aComponent)
{
    component = &aComponent;
    cur_S=0;
    tot_S=component->GetLength();
}

inline double ComponentIntegrator::GetIntegratedLength () const
{
    return cur_S;
}

inline double ComponentIntegrator::GetRemainingLength () const
{
    double s = tot_S-cur_S;
    return fequal(s,0) ? 0 : s;
}

inline bool ComponentIntegrator::IsValidStep (double ds) const
{
    return GetValidStepRange()(ds);
}

inline double ComponentIntegrator::Track(double ds)
{
    if(AtEntrance())
        TrackEntrance();

    TrackStep(ds);
    cur_S+=ds;

    if(AtExit())
        TrackExit();

    return GetRemainingLength();
}

inline bool ComponentIntegrator::AtEntrance() const
{
    return cur_S==0;
}

inline bool ComponentIntegrator::AtExit(double s) const
{
    return fequal(s+cur_S,tot_S);
}

inline FloatRange ComponentIntegrator::GetValidStepRange() const
{
    return FloatRange(0,GetRemainingLength());
}

#endif
