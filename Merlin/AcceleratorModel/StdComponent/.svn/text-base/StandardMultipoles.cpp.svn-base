//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\AcceleratorModel\StdComponent\StandardMultipoles.cpp
 * last modified 16/05/02 11:10:34
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



// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// StandardMultipoles
#include "AcceleratorModel/StdComponent/StandardMultipoles.h"


#define _RMC1(n)  RectMultipole(id,len,n,dnB)
#define _RMC2(n)  RectMultipole(id,len,n,B,r0)
#define _RMSC1(n)  RectMultipole(id,len,n,dnB,true)
#define _RMSC2(n)  RectMultipole(id,len,n,B,r0,true)
#define _MPT _PREPTRACK(aTracker,RectMultipole);
#define _RID return ID;
#define _CP(type) return new type (*this);


// Class Quadrupole

const int Quadrupole::ID = UniqueIndex();

Quadrupole::Quadrupole (const string& id, double len, double dnB)
        : _RMC1(1)
{
}

Quadrupole::Quadrupole (const string& id, double len, double B, double r0)
        : _RMC2(1)
{
}



void Quadrupole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int Quadrupole::GetIndex () const
{
    _RID
}

const string& Quadrupole::GetType () const
{
    _TYPESTR(Quadrupole)
}

ModelElement* Quadrupole::Copy () const
{
    _CP(Quadrupole)
}

// Class Sextupole

const int Sextupole::ID = UniqueIndex();

Sextupole::Sextupole (const string& id, double len, double dnB)
        : _RMC1(2)
{
}

Sextupole::Sextupole (const string& id, double len, double B, double r0)
        : _RMC2(2)
{
}



void Sextupole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int Sextupole::GetIndex () const
{
    _RID
}

const string& Sextupole::GetType () const
{
    _TYPESTR(Sextupole);
}

ModelElement* Sextupole::Copy () const
{
    _CP(Sextupole)
}

// Class SkewQuadrupole

const int SkewQuadrupole::ID = UniqueIndex();

SkewQuadrupole::SkewQuadrupole (const string& id, double len, double dnB)
        : _RMSC1(1)
{
}

SkewQuadrupole::SkewQuadrupole (const string& id, double len, double B, double r0)
        : _RMSC2(1)
{
}



void SkewQuadrupole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int SkewQuadrupole::GetIndex () const
{
    _RID
}

const string& SkewQuadrupole::GetType () const
{
    _TYPESTR(SkewQuadrupole)
}

ModelElement* SkewQuadrupole::Copy () const
{
    _CP(SkewQuadrupole)
}

// Class Octupole

const int Octupole::ID = UniqueIndex();

Octupole::Octupole (const string& id, double len, double dnB)
        : _RMC1(3)
{
}

Octupole::Octupole (const string& id, double len, double B, double r0)
        : _RMC2(3)
{
}



void Octupole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int Octupole::GetIndex () const
{
    _RID
}

const string& Octupole::GetType () const
{
    _TYPESTR(Octupole)
}

ModelElement* Octupole::Copy () const
{
    _CP(Octupole)
}

// Class SkewSextupole

const int SkewSextupole::ID = UniqueIndex();

SkewSextupole::SkewSextupole (const string& id, double len, double dnB)
        : _RMSC1(2)
{
}

SkewSextupole::SkewSextupole (const string& id, double len, double B, double r0)
        : _RMSC2(2)
{
}



void SkewSextupole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int SkewSextupole::GetIndex () const
{
    _RID
}

const string& SkewSextupole::GetType () const
{
    _TYPESTR(SkewSextupole)
}

ModelElement* SkewSextupole::Copy () const
{
    _CP(SkewSextupole)
}

// Class Decapole
const int Decapole::ID = UniqueIndex();

Decapole::Decapole (const string& id, double len, double dnB)
        : _RMC1(4)
{
}

Decapole::Decapole (const string& id, double len, double B, double r0)
        : _RMC2(4)
{
}



void Decapole::PrepareTracker (ComponentTracker& aTracker)
{
    _MPT
}

int Decapole::GetIndex () const
{
    _RID
}

const string& Decapole::GetType () const
{
    _TYPESTR(Decapole)
}

ModelElement* Decapole::Copy () const
{
    _CP(Decapole)
}

