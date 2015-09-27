/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		21.09.15 Haroon Rafique
// Modified:		
// Last Edited: 21.09.15 Haroon Rafique
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/StdComponent/CrabMarker.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int CrabMarker::ID = UniqueIndex();

CrabMarker::CrabMarker (const string& id, double len)
		: Drift(id,len)
{}

CrabMarker::CrabMarker (const string& id, double len, double mux, double muy)
		: Drift(id,len)
{
		SetMuX(mux);
		SetMuY(muy);
}

const string& CrabMarker::GetType () const
{
    _TYPESTR(CrabMarker);
}

int CrabMarker::GetIndex () const
{
    return  ID;
}

void CrabMarker::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Drift);
}

void CrabMarker::RotateY180 ()
{
    // nothing to do
}

ModelElement* CrabMarker::Copy () const
{
    return new CrabMarker(*this);
}

