/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "EnergyAdjustmentPolicy.h"

void EnergyAdjustmentPolicy::SetKlystrons(const KlystronArray& klys)
{
	theKlystrons = klys;
	defkvals.clear();
	for(size_t k = 0; k < theKlystrons.size(); k++)
	{
		defkvals.push_back(theKlystrons[k]->GetVoltagePhasor());
	}
}

void EnergyAdjustmentPolicy::RestoreKlystrons()
{
	for(size_t k = 0; k < theKlystrons.size(); k++)
	{
		theKlystrons[k]->SetVoltagePhasor(defkvals[k]);
	}
}
