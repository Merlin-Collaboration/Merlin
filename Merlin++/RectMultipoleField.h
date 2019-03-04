/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RectMultipoleField_h
#define RectMultipoleField_h 1

#include "merlin_config.h"
#include "TemplateComponents.h"
#include "MultipoleField.h"
#include "RectangularGeometry.h"

/**
 *	A magnetic multipole field referenced to a rectangular
 *	geometry.
 */

typedef TAccCompGF<RectangularGeometry, MultipoleField> RectMultipoleField;

#endif
