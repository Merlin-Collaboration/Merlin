/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ComponentTracker.h"
#include "Marker.h"

const int Marker::ID = UniqueIndex();

void Marker::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent);
}

int Marker::GetIndex() const
{
	return ID;
}

const string& Marker::GetType() const
{
	_TYPESTR(Marker)
}

//FIXME
ModelElement* Marker::Copy() const
{
	return new Marker(*this);
}

void Marker::RotateY180()
{
	// Nothing to do
}
