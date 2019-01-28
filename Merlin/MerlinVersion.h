/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <string>

/** \file MerlinVersion.h
 * Get information about the Merlin version
 */

/// Get the full Merlin version number
std::string merlin_version();

/// Get the Merlin version major number
int merlin_version_major();

/// Get the Merlin version minor number
int merlin_version_minor();

/// Get the Merlin version patch number
int merlin_version_patch();

/// Get the Merlin git tag (empty string if not compiled from git)
std::string merlin_version_git();

/// Get the Merlin version as printable text
std::string merlin_version_info();
