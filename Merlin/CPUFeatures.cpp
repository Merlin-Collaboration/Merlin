/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <string>
#include "CPUFeatures.h"
#include <cmath>

namespace CPUFeatures
{

void CheckCPUFeatures();
unsigned int GetCPUFeatures1();
unsigned int GetCPUFeatures2();
std::string GetCPUName();

void CheckCPUFeatures()
{
#ifdef __x86_64__
	std::string CPUName = GetCPUName();

	unsigned int ecxf = GetCPUFeatures1();
	unsigned int edxf = GetCPUFeatures2();

	std::cout << "Running on: " << CPUName << std::endl;

	//bool GotSSE2 = false;
	//bool GotAVX = false;

	//SSE2
	std::cout << "SSE2:\t";
	if((edxf >> 26) & 0x1)
	{
		std::cout << "Supported" << std::endl;
		//GotSSE2 = true;
	}
	else
	{
		std::cout << "Not suported" << std::endl;
	}

	//AVX
	std::cout << "AVX:\t";
	if((ecxf >> 28) & 0x1)
	{
		std::cout << "Supported" << std::endl;
		//GotAVX = true;
	}
	else
	{
		std::cout << "Not suported" << std::endl;
	}

#else
	std::cerr << "Currently only supported on X86_64" << std::endl;
#endif
}

/*
 * Get info from the ecx register
 */
unsigned int GetCPUFeatures1()
{
#ifdef __x86_64__
	unsigned int ecx;

	asm ("cpuid"
	: "=c" (ecx)
	: "a" (1)
	: "%ebx", "edx"
	);

	return ecx;
#else
	std::cerr << "Currently only supported on X86_64" << std::endl;
	return 0;
#endif
}

/*
 * Get info from the edx register
 */
unsigned int GetCPUFeatures2()
{
#ifdef __x86_64__
	unsigned int edx;

	asm ("cpuid"
	: "=d" (edx)
	: "a" (1)
	: "%ebx", "ecx"
	);

	return edx;
#else
	std::cerr << "Currently only supported on X86_64" << std::endl;
	return 0;
#endif
}

/*
 * Gets the CPU Name string
 * 3 calls are needed with eax = 0x80000002, 0x80000003 and 0x80000004
 */
std::string GetCPUName()
{
#ifdef __x86_64__
	std::string CPUNameString;
	char eax[4], ebx[4], ecx[4], edx[4];

	for(int j = 0; j < 3; j++)
	{
		asm ("cpuid"
		: "=a" (eax),
		"=b" (ebx),
		"=c" (ecx),
		"=d" (edx)
		: "a" (0x80000002 + j)
		);

		for(int i = 0; i < 4; i++)
		{
			CPUNameString.push_back(eax[i]);
		}
		for(int i = 0; i < 4; i++)
		{
			CPUNameString.push_back(ebx[i]);
		}
		for(int i = 0; i < 4; i++)
		{
			CPUNameString.push_back(ecx[i]);
		}
		for(int i = 0; i < 4; i++)
		{
			CPUNameString.push_back(edx[i]);
		}
	}

	return CPUNameString;
#else
	std::cerr << "Currently only supported on X86_64" << std::endl;
	return "UNKNOWN CPU";
#endif
}

#ifdef LIBNUMA
bool CheckNUMA()
{
	if(numa_available() < 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}

void PrintNUMAInfo()
{
	if(CheckNUMA())
	{
		std::cout << "NUMA info:" << std::endl;
		std::cout << "The kernel supports up to " << numa_num_possible_nodes() << " nodes." << std::endl;

		int n_nodes = numa_max_node();

		std::cout << "We have " << n_nodes + 1 << " NUMA nodes." << std::endl;

		std::cout << "We have " << numa_num_configured_cpus() << " CPU threads in total." << std::endl;

		std::cout << "We can use " << numa_num_task_cpus() << " CPU threads." << std::endl;
		std::cout << "We can allocate memory on " << numa_num_task_nodes() << " nodes." << std::endl;

		for(int n = 0; n < n_nodes + 1; n++)
		{
			std::cout << "There is " << numa_node_size(n, NULL) / pow(2, 20) << "Mb in total on node " << n
					  << std::endl;
		}
	}
	else
	{
		std::cout << "NUMA failure" << std::endl;
	}
}

#endif

} //End namespace
