/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ComponentTracker.h"
#include "StandardMultipoles.h"

#define _RMC1(n)  RectMultipole(id, len, n, dnB)
#define _RMC2(n)  RectMultipole(id, len, n, B, r0)
#define _RMSC1(n)  RectMultipole(id, len, n, dnB, true)
#define _RMSC2(n)  RectMultipole(id, len, n, B, r0, true)
#define _MPT _PREPTRACK(aTracker, RectMultipole);
#define _RID return ID;
#define _CP(type) return new type(*this);

// Class Quadrupole
const int Quadrupole::ID = UniqueIndex();

Quadrupole::Quadrupole(const string& id, double len, double dnB) :
	_RMC1(1)
{
}

Quadrupole::Quadrupole(const string& id, double len, double B, double r0) :
	_RMC2(1)
{
}

void Quadrupole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int Quadrupole::GetIndex() const
{
	_RID
}

const string& Quadrupole::GetType() const
{
	_TYPESTR(Quadrupole)
}

ModelElement* Quadrupole::Copy() const
{
	_CP(Quadrupole)
}

// Class Sextupole

const int Sextupole::ID = UniqueIndex();

Sextupole::Sextupole(const string& id, double len, double dnB) :
	_RMC1(2)
{
}

Sextupole::Sextupole(const string& id, double len, double B, double r0) :
	_RMC2(2)
{
}

void Sextupole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int Sextupole::GetIndex() const
{
	_RID
}

const string& Sextupole::GetType() const
{
	_TYPESTR(Sextupole);
}

ModelElement* Sextupole::Copy() const
{
	_CP(Sextupole)
}

// Class SkewQuadrupole

const int SkewQuadrupole::ID = UniqueIndex();

SkewQuadrupole::SkewQuadrupole(const string& id, double len, double dnB) :
	_RMSC1(1)
{
}

SkewQuadrupole::SkewQuadrupole(const string& id, double len, double B, double r0) :
	_RMSC2(1)
{
}

void SkewQuadrupole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int SkewQuadrupole::GetIndex() const
{
	_RID
}

const string& SkewQuadrupole::GetType() const
{
	_TYPESTR(SkewQuadrupole)
}

ModelElement* SkewQuadrupole::Copy() const
{
	_CP(SkewQuadrupole)
}

// Class Octupole

const int Octupole::ID = UniqueIndex();

Octupole::Octupole(const string& id, double len, double dnB) :
	_RMC1(3)
{
}

Octupole::Octupole(const string& id, double len, double B, double r0) :
	_RMC2(3)
{
}

void Octupole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int Octupole::GetIndex() const
{
	_RID
}

const string& Octupole::GetType() const
{
	_TYPESTR(Octupole)
}

ModelElement* Octupole::Copy() const
{
	_CP(Octupole)
}

// Class SkewSextupole

const int SkewSextupole::ID = UniqueIndex();

SkewSextupole::SkewSextupole(const string& id, double len, double dnB) :
	_RMSC1(2)
{
}

SkewSextupole::SkewSextupole(const string& id, double len, double B, double r0) :
	_RMSC2(2)
{
}

void SkewSextupole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int SkewSextupole::GetIndex() const
{
	_RID
}

const string& SkewSextupole::GetType() const
{
	_TYPESTR(SkewSextupole)
}

ModelElement* SkewSextupole::Copy() const
{
	_CP(SkewSextupole)
}

// Class Decapole
const int Decapole::ID = UniqueIndex();

Decapole::Decapole(const string& id, double len, double dnB) :
	_RMC1(4)
{
}

Decapole::Decapole(const string& id, double len, double B, double r0) :
	_RMC2(4)
{
}

void Decapole::PrepareTracker(ComponentTracker& aTracker)
{
	_MPT
}

int Decapole::GetIndex() const
{
	_RID
}

const string& Decapole::GetType() const
{
	_TYPESTR(Decapole)
}

ModelElement* Decapole::Copy() const
{
	_CP(Decapole)
}
