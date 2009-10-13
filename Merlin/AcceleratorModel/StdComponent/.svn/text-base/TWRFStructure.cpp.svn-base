/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

// TWRFStructure
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// TWRFfield
#include "AcceleratorModel/StdField/TWRFfield.h"

const int TWRFStructure::ID = UniqueIndex();

TWRFStructure::TWRFStructure (const string& id, double len, double f, double Epk, double phi)
: RFStructure(id,len,new TWRFfield(f,Epk,phi))
{}

TWRFStructure::TWRFStructure (const TWRFStructure& rhs)
:RFStructure(rhs.GetName(),rhs.GetLength(),new TWRFfield(static_cast<const TWRFfield&>(rhs.GetField())))
{}

const string& TWRFStructure::GetType () const
{
	_TYPESTR(TWRFStructure)
}

int TWRFStructure::GetIndex () const
{
	return ID;
}

void TWRFStructure::PrepareTracker (ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker,AcceleratorComponent)
}

void TWRFStructure::RotateY180 ()
{
	double E = GetField().GetAmplitude();
	GetField().SetAmplitude(-E);
}

ModelElement* TWRFStructure::Copy () const
{
	return new TWRFStructure(*this);
}

