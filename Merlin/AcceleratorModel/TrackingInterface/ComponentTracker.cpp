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

#include "merlin_config.h"
#include <cassert>
#include <algorithm>
#include <map>
#include <limits>
#include <iostream>
#include "stdext/deleters.h"
#include "NumericalUtils/utils.h"
// ComponentIntegrator
#include "AcceleratorModel/TrackingInterface/ComponentIntegrator.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// Marker
#include "AcceleratorModel/StdComponent/Marker.h"

//	The onlt default integrator supplied is that for a
//	Marker element, which are zero-length pseudo-components.
//	The default integrator does nothing.

class DefaultMarkerIntegrator : public ComponentIntegrator
{
public:

    //	Performs no action.
    virtual void TrackStep (double ds);

    //	Returns the component index for this integrator.
    virtual int GetComponentIndex () const;
};

ComponentTracker::IntegratorSet::~IntegratorSet ()
{
    //	std::for_each(itsMap.begin(),itsMap.end(),map_deleter<IndexType,Integrator>());
    for(IMap::iterator i=itsMap.begin();i!=itsMap.end();i++) {
        delete (*i).second;
    }
}

bool ComponentTracker::IntegratorSet::Add (ComponentIntegrator* intg)
{
    using namespace std;

    pair<IMap::iterator,bool> rt=itsMap.insert(
                                     IMap::value_type(intg->GetComponentIndex(),intg));

    if(!rt.second) { // override
        delete (*rt.first).second;
        (*rt.first).second = intg;
    }
    return !rt.second;
}

ComponentIntegrator* ComponentTracker::IntegratorSet::GetIntegrator (int n)
{
    IMap::iterator i = itsMap.find(n);
    return i==itsMap.end() ? 0 : (*i).second;
}

void DefaultMarkerIntegrator::TrackStep (double ds)
{
    assert(ds==0);
}

int DefaultMarkerIntegrator::GetComponentIndex () const
{
    return Marker::ID;
}

ComponentTracker::ComponentTracker ()
        : itsState(undefined),iSet(new IntegratorSet)
{
    Register(new DefaultMarkerIntegrator());
}

ComponentTracker::ComponentTracker (IntegratorSet* anIS)
        : itsState(undefined),iSet(anIS)
{}

ComponentTracker::~ComponentTracker ()
{
    if(iSet)
        delete iSet;
}

void ComponentTracker::Track ()
{
    assert(itsState==initialised);
    integrator->TrackAll();
    itsState=finished;
}

double ComponentTracker::TrackStep (double ds)
{
    assert(itsState==initialised || itsState==tracking);
    assert(integrator->IsValidStep(ds));

    double sToExit=integrator->Track(ds);
    itsState = fequal(sToExit,0) ? finished : tracking;
    return itsState;
}

void ComponentTracker::Reset ()
{
    integrator = 0;
    itsState = undefined;
}

double ComponentTracker::GetRemainingLength () const
{
    assert(itsState!=undefined);
    return integrator->GetRemainingLength();
}

double ComponentTracker::GetIntegratedLength () const
{
    assert(itsState!=undefined);
    return integrator->GetIntegratedLength();
}

bool ComponentTracker::SelectIntegrator (int index, AcceleratorComponent& component)
{
    assert((itsState==undefined)||(itsState==finished));

    integrator = iSet->GetIntegrator(index);
    if(!integrator)
        return false;
    integrator->SetCurrentComponent(component);
    InitialiseIntegrator(integrator);
    return true;
}

void ComponentTracker::InitialiseIntegrator (ComponentIntegrator*)
{
    itsState = initialised;
}

bool ComponentTracker::Register (ComponentIntegrator* intg)
{
    if(iSet==0)
        iSet = new IntegratorSet();

    return iSet->Add(intg);
}




