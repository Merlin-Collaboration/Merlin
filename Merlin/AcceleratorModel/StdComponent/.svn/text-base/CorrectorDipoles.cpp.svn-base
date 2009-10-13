/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

// CorrectorDipoles
#include "AcceleratorModel/StdComponent/CorrectorDipoles.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int XCor::ID = UniqueIndex();
const int YCor::ID = UniqueIndex();

ModelElement* XCor::Copy () const
{
    return new XCor(*this);
}

int XCor::GetIndex () const
{
    return ID;
}

const string& XCor::GetType () const
{
    _TYPESTR(XCor)
}

void XCor::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,RectMultipole)
}

ModelElement* YCor::Copy () const
{
    return new YCor(*this);
}

int YCor::GetIndex () const
{
    return ID;
}

const string& YCor::GetType () const
{
    _TYPESTR(YCor)
}

void YCor::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,RectMultipole)
}

