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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/AcceleratorGeometry.h"

Transform3D AcceleratorGeometry::GetTotalGeometryTransform () const
{
    Extent extnt = GetGeometryExtent();
    return GetGeometryTransform(extnt.first,extnt.second);
}

