/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RotationType_h
#define RotationType_h 1

#include "merlin_config.h"

/**
 *	Type of rotation. Can be ident, xrot, yrot, zrot or mixed
 */

typedef enum
{
	ident,
	xrot,
	yrot,
	zrot,
	mixed

} RotationType;

#endif
