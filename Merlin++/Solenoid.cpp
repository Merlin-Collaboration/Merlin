/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ComponentTracker.h"

#include "Solenoid.h"

// Class Solenoid
const int Solenoid::ID = UniqueIndex();

Solenoid::Solenoid(const std::string& id, double len, double Bz) :
	SimpleSolenoid(id, new RectangularGeometry(len), new BzField(Bz))
{
}

void Solenoid::RotateY180()
{
	BzField& field = GetField();
	field.SetStrength(-field.GetStrength());
}

const string& Solenoid::GetType() const
{
	_TYPESTR(Solenoid);
}

ModelElement* Solenoid::Copy() const
{
	return new Solenoid(*this);
}

int Solenoid::GetIndex() const
{
	return Solenoid::ID;
}

void Solenoid::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent);
}
