//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\AcceleratorModel\StdComponent\Solenoid.cpp
 * last modified 10/12/01 16:41:41
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


#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

// Solenoid
#include "AcceleratorModel/StdComponent/Solenoid.h"


// Class Solenoid

const int Solenoid::ID = UniqueIndex();

Solenoid::Solenoid (const std::string& id, double len, double Bz)
        : SimpleSolenoid(id,new RectangularGeometry(len),new BzField(Bz))
{
}



void Solenoid::RotateY180 ()
{
    BzField& field = GetField();
    field.SetStrength(-field.GetStrength());
}

const string& Solenoid::GetType () const
{
    _TYPESTR(Solenoid);
}

ModelElement* Solenoid::Copy () const
{
    return new Solenoid(*this);
}

int Solenoid::GetIndex () const
{
    return Solenoid::ID;
}

void Solenoid::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

