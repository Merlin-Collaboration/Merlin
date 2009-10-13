/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// EMField
#include "AcceleratorModel/EMField.h"
// AcceleratorGeometry
#include "AcceleratorModel/AcceleratorGeometry.h"
// Aperture
#include "AcceleratorModel/Aperture.h"
// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"

const int AcceleratorComponent::ID = UniqueIndex();

AcceleratorComponent::~AcceleratorComponent ()
{
    if(itsGeometry) delete itsGeometry;
    if(itsField) delete itsField;
}

int AcceleratorComponent::GetIndex () const
{
    return ID;
}

double AcceleratorComponent::GetLength () const
{
    return itsGeometry ? itsGeometry->GetGeometryLength() : 0;
}

void AcceleratorComponent::PrepareTracker (ComponentTracker& aTracker)
{
    if(!aTracker.SelectIntegrator(AcceleratorComponent::ID,*this))
        throw ComponentTracker::UnknownComponent();
}

int AcceleratorComponent::UniqueIndex ()
{
    static int ID_count = 0;
    return ID_count++;
}

