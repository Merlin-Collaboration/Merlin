/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TWRFStructure.h"
#include "ComponentTracker.h"
#include "TWRFfield.h"

const int TWRFStructure::ID = UniqueIndex();

TWRFStructure::TWRFStructure(const string& id, double len, double f, double Epk, double phi) :
	RFStructure(id, len, new TWRFfield(f, Epk, phi))
{
}

TWRFStructure::TWRFStructure(const TWRFStructure& rhs) :
	RFStructure(rhs.GetName(), rhs.GetLength(), new TWRFfield(static_cast<const TWRFfield&>(rhs.GetField())))
{
}

const string& TWRFStructure::GetType() const
{
	_TYPESTR(TWRFStructure)
}

int TWRFStructure::GetIndex() const
{
	return ID;
}

void TWRFStructure::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent)
}

void TWRFStructure::RotateY180()
{
	double E = GetField().GetAmplitude();
	GetField().SetAmplitude(-E);
}

ModelElement* TWRFStructure::Copy() const
{
	return new TWRFStructure(*this);
}
