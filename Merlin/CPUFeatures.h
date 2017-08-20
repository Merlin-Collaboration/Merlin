#ifndef _CPUFeatures_h_
#define _CPUFeatures_h_ 1

#include <iostream>
#include <string>

#ifdef LIBNUMA
#include <numa.h>
#endif

namespace CPUFeatures
{

void CheckCPUFeatures();
unsigned int GetCPUFeatures1();
unsigned int GetCPUFeatures2();
std::string GetCPUName();

/**
* NUMA features (Non Uniform Memory Access)
*
* @see `man numa` on linux for details
*/
#ifdef LIBNUMA
bool CheckNUMA();
void PrintNUMAInfo();
#endif

}//End namespace

#endif

