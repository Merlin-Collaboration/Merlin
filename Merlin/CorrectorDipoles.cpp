/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "CorrectorDipoles.h"
#include "ComponentTracker.h"

const int XCor::ID = UniqueIndex();
const int YCor::ID = UniqueIndex();

ModelElement* XCor::Copy() const
{
	return new XCor(*this);
}

int XCor::GetIndex() const
{
	return ID;
}

const string& XCor::GetType() const
{
	_TYPESTR(XCor)
}

void XCor::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, RectMultipole)
}

ModelElement* YCor::Copy() const
{
	return new YCor(*this);
}

int YCor::GetIndex() const
{
	return ID;
}

const string& YCor::GetType() const
{
	_TYPESTR(YCor)
}

void YCor::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, RectMultipole)
}
