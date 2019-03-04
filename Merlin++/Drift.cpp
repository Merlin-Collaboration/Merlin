/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Drift.h"
#include "ComponentTracker.h"

const int Drift::ID = UniqueIndex();

Drift::Drift(const string& id, double len) :
	TAccCompG<RectangularGeometry>(id, new RectangularGeometry(len))
{
}

const string& Drift::GetType() const
{
	_TYPESTR(Drift);
}

int Drift::GetIndex() const
{
	return ID;
}

void Drift::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent);
}

void Drift::RotateY180()
{
	// nothing to do
}

ModelElement* Drift::Copy() const
{
	return new Drift(*this);
}
