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

#ifndef RectMultipoleField_h
#define RectMultipoleField_h 1

#include "merlin_config.h"
// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"
// RectangularGeometry
#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"

//	A magnetic multipole field referenced to a rectangular
//	geometry.

typedef TAccCompGF< RectangularGeometry,MultipoleField > RectMultipoleField;

#endif
