
/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (2000)
* 
* file Merlin\AcceleratorModel\StdComponent\Spoiler.h
* last modified 04/04/01 15:25:43
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
*
* Copyright (c) 2000 by The Merlin Collaboration.  
* ALL RIGHTS RESERVED. 
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

#ifndef Spoiler_h
#define Spoiler_h 1

#include "merlin_config.h"

#include "AcceleratorModel/StdComponent/Drift.h"

class ComponentTracker;

// A spoiler represents a scattering element in the beamline. Spoiler objects
// are optically drifts, but when associated with an Aperture, they have the
// material property of radiation length, which is used by specific scattering
// process to approximate the interaction of particles with the spoiler material.

class Spoiler : public Drift
{
public:

    Spoiler (const string& id, double len, double radLength);

    // Returns the length of the spoiler in units of its
    // radiation length
    double GetNumRadLengths() const {
        return GetLength()/Xr;
    }
	
	// Returns the material radiation length (meter)
	double GetMaterialRadiationLength() const {
		return Xr;
	}

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Unique index for an Accelerator component.
    static const int ID;

private:

    double Xr;
};

#endif
