/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Bunch.h"
#include "ComponentTracker.h"
#include "Monitor.h"

bool Monitor::all_inactive = false;
const int Monitor::ID = UniqueIndex();

Monitor::Monitor(const string& id, double len, double mpt) :
	TAccCompG<RectangularGeometry>(id, new RectangularGeometry(len)), mp(mpt), active(true)
{
}

Monitor::~Monitor()
{
}

void Monitor::MakeMeasurement(const Bunch&)
{
}

void Monitor::SetMeasurementPt(double mpt)
{
	//GetGeometry().CheckBounds(mpt); // might throw
	mp = mpt;
}

double Monitor::GetMeasurementPt() const
{
	return mp;
}

void Monitor::RotateY180()
{
	reflected = !reflected;
}

int Monitor::GetIndex() const
{
	return ID;
}

const string& Monitor::GetType() const
{
	_TYPESTR(Monitor)
}

void Monitor::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent)
}

ModelElement* Monitor::Copy() const
{
	return new Monitor(*this);
}
