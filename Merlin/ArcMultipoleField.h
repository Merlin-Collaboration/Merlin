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

#ifndef ArcMultipoleField_h
#define ArcMultipoleField_h 1

#include "merlin_config.h"
// TemplateComponents
#include "TemplateComponents.h"
// MultipoleField
#include "MultipoleField.h"
// ArcGeometry
#include "ArcGeometry.h"

/**
*	A magnetic multipole field referenced to an arc geometry.
*/
typedef TAccCompGF< ArcGeometry,MultipoleField > ArcMultipoleField;

#endif
