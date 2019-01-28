/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TransverseRFStructure.h"
#include "ComponentTracker.h"

const int TransverseRFStructure::ID = UniqueIndex();

TransverseRFStructure::TransverseRFStructure(const string& id, double len, double f, double Epk, double phi, double
	theta) :
	RFStructure(id, len, new TransverseRFfield(f, Epk, phi, theta))
{
}

TransverseRFStructure::TransverseRFStructure(const TransverseRFStructure& rhs) :
	RFStructure(rhs.GetName(), rhs.GetLength(), new TransverseRFfield(static_cast<const
		TransverseRFfield&>(rhs.GetField())))
{
}

const string& TransverseRFStructure::GetType() const
{
	_TYPESTR(TransverseRFStructure)
}

int TransverseRFStructure::GetIndex() const
{
	return ID;
}

void TransverseRFStructure::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent)
}

void TransverseRFStructure::RotateY180()
{
	double E = GetField().GetAmplitude();
	double t = GetFieldOrientation();
	GetField().SetAmplitude(-E);
	SetFieldOrientation(-t);
}

ModelElement* TransverseRFStructure::Copy() const
{
	return new TransverseRFStructure(*this);
}
