/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "RectMultipole.h"
#include "ComponentTracker.h"

const int RectMultipole::ID = UniqueIndex();

RectMultipole::RectMultipole(const string& id, double length, int npole, double b, double r0, bool skew) :
	TAccCompGF<RectangularGeometry, MultipoleField>(id, new RectangularGeometry(length),
		new MultipoleField(npole, b, r0, skew))
{
}

RectMultipole::RectMultipole(const string& id, double len, int np, double b, bool skew) :
	TAccCompGF<RectangularGeometry, MultipoleField>(id, new RectangularGeometry(len),
		new MultipoleField(np, b, skew))
{
}

int RectMultipole::GetIndex() const
{
	return ID;
}

void RectMultipole::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent);
}

void RectMultipole::RotateY180()
{
	GetField().RotateY180();
}
