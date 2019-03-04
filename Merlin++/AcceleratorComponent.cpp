/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ComponentTracker.h"
#include "EMField.h"
#include "AcceleratorGeometry.h"
#include "Aperture.h"
#include "AcceleratorComponent.h"

const int AcceleratorComponent::ID = UniqueIndex();

AcceleratorComponent::~AcceleratorComponent()
{
	if(itsGeometry)
	{
		delete itsGeometry;
	}

	if(itsField)
	{
		delete itsField;
	}
}

int AcceleratorComponent::GetIndex() const
{
	return ID;
}

double AcceleratorComponent::GetLength() const
{
	return itsGeometry ? itsGeometry->GetGeometryLength() : 0;
}

void AcceleratorComponent::PrepareTracker(ComponentTracker& aTracker)
{
	if(!aTracker.SelectIntegrator(AcceleratorComponent::ID, *this))
	{
		throw ComponentTracker::UnknownComponent();
	}
}

int AcceleratorComponent::UniqueIndex()
{
	static int ID_count = 0;
	return ID_count++;
}
