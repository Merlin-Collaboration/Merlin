//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (1999)
* 
* file Merlin\BeamDynamics\ParticleTracking\ParticleMapComponent.h
* last modified 26/09/02 15:18:13
*/

/*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
* Copyright (c) 2001 by The Merlin Collaboration.
* - ALL RIGHTS RESERVED - 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/


#ifndef ParticleMapComponent_h
#define ParticleMapComponent_h 1

#include "merlin_config.h"


// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"
// ParticleMap
#include "BeamDynamics/ParticleTracking/ParticleMap.h"

namespace ParticleTracking {

class ParticleBunch;






//	A special AcceleratorComponent that allows an arbitrary
//	map (ParticleMap) to be placed into the accelerator
//	model.










class ParticleMapComponent : public AcceleratorComponent
{
public:


    ParticleMapComponent (const std::string& id, ParticleMap* pmap, double intB2ds = 0);




    //	Return the type string for the element.
    virtual const string& GetType () const;


    //	Virtual constructor.
    virtual ModelElement* Copy () const;


    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;


    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();


    ParticleBunch& Apply (ParticleBunch& bunch) const;


    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);


    //	Returns the integral of B^2 for synchrotron radiation
    //	applications
    double GetIntB2ds () const;

    // Data Members for Class Attributes


    //	Unique index for an Accelerator component.

    static const int ID;


protected:
private:
    // Data Members for Class Attributes



    double ib2;


    // Data Members for Associations




    ParticleMap* itsMap;


private:
};

// Class ParticleMapComponent




inline ParticleBunch& ParticleMapComponent::Apply (ParticleBunch& bunch) const
{

    return itsMap->Apply(bunch);

}


inline double ParticleMapComponent::GetIntB2ds () const
{

    return ib2;

}




}; //end namespace ParticleTracking

#endif
