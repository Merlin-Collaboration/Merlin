#include "utility/MerlinProfile.h"

#include <map>
#include <string>
std::map<std::string, timespec> MerlinProfile::pData;
std::map<std::string, timespec> MerlinProfile::StartTime;

// To test whether MERLIN_PROFILE was enabled when libmerlin was built
#ifdef MERLIN_PROFILE
bool MerlinProfile::IsEnabled(){return true;}
#else
bool MerlinProfile::IsEnabled(){return false;}
#endif
