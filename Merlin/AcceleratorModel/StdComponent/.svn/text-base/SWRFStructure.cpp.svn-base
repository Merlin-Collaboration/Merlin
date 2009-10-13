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

#include "AcceleratorModel/StdComponent/SWRFStructure.h"
#include "AcceleratorModel/StdField/SWRFfield.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

namespace {
inline double Wavelength(double f) {
    using PhysicalConstants::SpeedOfLight;
    return SpeedOfLight/f;
}
};

// Class SWRFStructure

const int SWRFStructure::ID = UniqueIndex();

SWRFStructure::SWRFStructure (const string& id, int ncells, double f, double E0, double phi)
: RFStructure(id,Wavelength(f)*ncells/2,new SWRFfield(f,E0,phi))
{}

SWRFStructure::SWRFStructure (const SWRFStructure& rhs)
: RFStructure(rhs.GetName(),rhs.GetLength(),new SWRFfield(static_cast<const SWRFfield&>(rhs.GetField())))
{}

const string& SWRFStructure::GetType () const
{
    _TYPESTR(SWRFStructure)
}

int SWRFStructure::GetIndex () const
{
    return ID;
}

void SWRFStructure::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent)
}

void SWRFStructure::RotateY180 ()
{
    double E = GetField().GetAmplitude();
    GetField().SetAmplitude(-E);
}

ModelElement* SWRFStructure::Copy () const
{
    return new SWRFStructure(*this);
}

