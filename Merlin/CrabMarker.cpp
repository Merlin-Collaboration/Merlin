/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "CrabMarker.h"
#include "ComponentTracker.h"

const int CrabMarker::ID = UniqueIndex();

CrabMarker::CrabMarker(const string& id, double len) :
	Drift(id, len)
{
}

CrabMarker::CrabMarker(const string& id, double len, double mux, double muy) :
	Drift(id, len)
{
	SetMuX(mux);
	SetMuY(muy);
}

const string& CrabMarker::GetType() const
{
	_TYPESTR(CrabMarker);
}

int CrabMarker::GetIndex() const
{
	return ID;
}

void CrabMarker::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, Drift);
}

void CrabMarker::RotateY180()
{
	// nothing to do
}

ModelElement* CrabMarker::Copy() const
{
	return new CrabMarker(*this);
}
