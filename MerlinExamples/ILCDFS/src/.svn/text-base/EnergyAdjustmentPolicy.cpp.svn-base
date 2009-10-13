/////////////////////////////////////////////////////////////////////////
// Abstract class EnergyAdjustmentPolicy implementation
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "EnergyAdjustmentPolicy.h"

void EnergyAdjustmentPolicy::SetKlystrons(const KlystronArray& klys)
{
	theKlystrons = klys;
	defkvals.clear();
	for(size_t k=0; k<theKlystrons.size(); k++)
		defkvals.push_back(theKlystrons[k]->GetVoltagePhasor());
}

void EnergyAdjustmentPolicy::RestoreKlystrons()
{
	for(size_t k=0; k<theKlystrons.size(); k++)
		theKlystrons[k]->SetVoltagePhasor(defkvals[k]);
}
