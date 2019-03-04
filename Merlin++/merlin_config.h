/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_merlin_config
#define _h_merlin_config
#include <string>
#include <cstring>

//Include for x86 vector functions
#ifdef __x86_64__
#include <immintrin.h>
#endif

//#if _MSC_VER >= 1400
//#define __TYPENAME__ typename
//#endif

#ifdef _MSC_VER
#pragma warning(disable: 4786)
#pragma warning(disable: 4660)
#pragma warning(disable: 4290)
#pragma warning(disable: 4018)
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#if _MSC_VER >= 1400
#define __TYPENAME__ typename
#else
#define __TYPENAME__
#endif

#else
#define __TYPENAME__ typename
#define _MAX(a, b) std::max(a, b)
#define _MIN(a, b) std::min(a, b)
#endif /* _MSC_VER */
//dk must be defined somewhere
#define _MAX(a, b) std::max(a, b)
#define _MIN(a, b) std::min(a, b)

//Define the floating point type used in Merlin
/** Use doubles (64-bit) */
typedef double Mfloat_t;

/** Use floats (32-bit) */
//typedef float Mfloat_t;

#endif /* _h_merlin_config */
