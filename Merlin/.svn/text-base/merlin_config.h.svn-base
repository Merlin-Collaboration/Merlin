// merlin_config.h
// configuration and platform/compiler dependent information
// Currently only needed for MSVC++

#ifndef _h_merlin_config
#define _h_merlin_config
#include <string>
#include <cstring>

#if _MSC_VER >= 1400
#define __TYPENAME__ typename
#endif

#ifdef _MSC_VER
#pragma warning(disable: 4786)
#pragma warning(disable: 4660)
#pragma warning(disable: 4290)
#pragma warning(disable: 4018)
#pragma warning(disable: 4244)
#pragma warning(disable: 4267)

#if _MSC_VER>=1400
#define __TYPENAME__ typename
#else
#define __TYPENAME__
#endif
#else
#define __TYPENAME__ typename
#define _MAX std::max
#define _MIN std::min
#endif /* _MSC_VER */
//dk must be defined somewhere
#define _MAX std::max
#define _MIN std::min

#endif /* _h_merlin_config */
