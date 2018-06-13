/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "SWRFStructure.h"
#include "SWRFfield.h"
#include "ComponentTracker.h"

namespace
{
inline double Wavelength(double f)
{
	using PhysicalConstants::SpeedOfLight;
	return SpeedOfLight / f;
}
}

// Class SWRFStructure

const int SWRFStructure::ID = UniqueIndex();

SWRFStructure::SWRFStructure(const string& id, int ncells, double f, double E0, double phi) :
	RFStructure(id, Wavelength(f) * ncells / 2, new SWRFfield(f, E0, phi))
{
}

SWRFStructure::SWRFStructure(const SWRFStructure& rhs) :
	RFStructure(rhs.GetName(), rhs.GetLength(), new SWRFfield(static_cast<const SWRFfield&>(rhs.GetField())))
{
}

const string& SWRFStructure::GetType() const
{
	_TYPESTR(SWRFStructure)
}

int SWRFStructure::GetIndex() const
{
	return ID;
}

void SWRFStructure::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent)
}

void SWRFStructure::RotateY180()
{
	double E = GetField().GetAmplitude();
	GetField().SetAmplitude(-E);
}

ModelElement* SWRFStructure::Copy() const
{
	return new SWRFStructure(*this);
}
