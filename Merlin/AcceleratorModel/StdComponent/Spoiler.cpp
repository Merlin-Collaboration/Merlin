
/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (2000)
* 
* file Merlin\AcceleratorModel\StdComponent\Spoiler.cpp
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


// Spoiler
#include "AcceleratorModel/StdComponent/Spoiler.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

// Class Spoiler
const int Spoiler::ID = UniqueIndex();

Spoiler::Spoiler (const string& id, double len, double radLength)
        : Drift(id,len),Xr(radLength)
{}

const string& Spoiler::GetType () const
{
    _TYPESTR(Spoiler);
}

int Spoiler::GetIndex () const
{
    return  ID;
}


void Spoiler::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Drift);
}

void Spoiler::RotateY180 ()
{
    // nothing to do
}


ModelElement* Spoiler::Copy () const
{
    return new Spoiler(*this);
}



