
/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
*
* Class library version 2.0 (2000)
*
* file Merlin\AcceleratorModel\StdComponent\Collimator.cpp
* last modified 04/04/01 15:25:43
*
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

#include "AcceleratorModel/StdComponent/Collimator.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

// Class Collimator
const int Collimator::ID = UniqueIndex();

Collimator::Collimator (const string& id, double len, double radLength)
        : Drift(id,len),Xr(radLength)
{scatter_at_this_collimator = true;}

Collimator::Collimator (const string& id, double len)
        : Drift(id,len),Xr(0)
{scatter_at_this_collimator = true;}

const string& Collimator::GetType () const
{
    _TYPESTR(Collimator);
}

int Collimator::GetIndex () const
{
    return  ID;
}

void Collimator::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Drift);
}

void Collimator::RotateY180 ()
{
    // nothing to do
}

ModelElement* Collimator::Copy () const
{
    return new Collimator(*this);
}



