/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		06.10.14 Haroon Rafique
// Modified:		
// Last Edited: 05.01.15 HR
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/StdComponent/HollowElectronLens.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int HollowElectronLens::ID = UniqueIndex();

HollowElectronLens::HollowElectronLens (const string& id, double len)
		: Drift(id,len)
{}

const string& HollowElectronLens::GetType () const
{
    _TYPESTR(HollowElectronLens);
}

int HollowElectronLens::GetIndex () const
{
    return  ID;
}


void HollowElectronLens::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Drift); //HR change to stop collimator tracking 29.10.13
    //_PREPTRACK(aTracker,AcceleratorComponent);
}

void HollowElectronLens::RotateY180 ()
{
    // nothing to do
}

void  HollowElectronLens::SetRmax (double rmax)
{
	Rmax = rmax;
}

void  HollowElectronLens::SetRmin (double rmin)
{
	Rmin = rmin;
}

ModelElement* HollowElectronLens::Copy () const
{
    return new HollowElectronLens(*this);
}

