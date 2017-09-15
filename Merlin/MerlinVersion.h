#include <string>

/** \file MerlinVersion.h
*Get information about the Merlin version
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
