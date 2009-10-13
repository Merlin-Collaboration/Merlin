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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////


#include "AcceleratorModel/StdComponent/TransverseRFStructure.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int TransverseRFStructure::ID = UniqueIndex();

TransverseRFStructure::TransverseRFStructure (const string& id, double len, double f, double Epk, double phi, double theta)
        : RFStructure(id,len,new TransverseRFfield(f,Epk,phi,theta))
{}

TransverseRFStructure::TransverseRFStructure (const TransverseRFStructure& rhs)
: RFStructure(rhs.GetName(),rhs.GetLength(),new TransverseRFfield(static_cast<const TransverseRFfield&>(rhs.GetField())))
{}

const string& TransverseRFStructure::GetType () const
{
    _TYPESTR(TransverseRFStructure)
}

int TransverseRFStructure::GetIndex () const
{
    return ID;
}

void TransverseRFStructure::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent)
}

void TransverseRFStructure::RotateY180 ()
{
    double E = GetField().GetAmplitude();
    double t = GetFieldOrientation();
    GetField().SetAmplitude(-E);
    SetFieldOrientation(-t);
}

ModelElement* TransverseRFStructure::Copy () const
{
    return new TransverseRFStructure(*this);
}

