/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Collimator.h"
#include "ComponentTracker.h"

// Class Collimator
const int Collimator::ID = UniqueIndex();

Collimator::Collimator(const string& id, double len) :
	Drift(id, len), Xr(0), scatter_at_this_collimator(true)
{

}

Collimator::Collimator(const string& id, double len, double radLength) :
	Drift(id, len), Xr(radLength), scatter_at_this_collimator(true)
{

}

Collimator::Collimator(const string& id, double len, Material* mat, double P0) :
	Drift(id, len), material(mat), scatter_at_this_collimator(true)
{

}

const string& Collimator::GetType() const
{
	_TYPESTR(Collimator);
}

int Collimator::GetIndex() const
{
	return ID;
}

void Collimator::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, Drift);
}

void Collimator::RotateY180()
{
	// nothing to do
}

ModelElement* Collimator::Copy() const
{
	return new Collimator(*this);
}
