/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _MerlinProfile_hpp_
#define _MerlinProfile_hpp_ 1

// Macros that do nothing if profiling not enabled
#ifdef MERLIN_PROFILE
#define MERLIN_PROFILE_ADD_PROCESS(s) MerlinProfile::AddProcess(s)
#define MERLIN_PROFILE_START_TIMER(s) MerlinProfile::StartProcessTimer(s)
#define MERLIN_PROFILE_END_TIMER(s) MerlinProfile::EndProcessTimer(s)
#else
#define MERLIN_PROFILE_ADD_PROCESS(s)
#define MERLIN_PROFILE_START_TIMER(s)
#define MERLIN_PROFILE_END_TIMER(s)
#endif

#include <string>
#include <ctime>
#include <iostream>
#include <map>

#ifdef MERLIN_PROFILE

#include <utility>
#include <algorithm>
#include <iomanip>

class MerlinProfile
{

//struct timespec {
//        time_t   tv_sec;        /* seconds */
//        long     tv_nsec;       /* nanoseconds */
//};

public:
//When adding a process want to make a new entry
	static void AddProcess(std::string ID)
	{
		std::map<std::string, timespec>::iterator Location = pData.begin();
		Location = pData.find(ID);
		if(Location == pData.end())
		{
			std::cout << "Adding Profile for: " << ID << std::endl;
			timespec pTime;
			pTime.tv_nsec = 0;
			pTime.tv_sec = 0;
			pData.insert(std::pair<std::string, timespec>(ID, pTime));
			StartTime.insert(std::pair<std::string, timespec>(ID, pTime));
		}
	}

//When removing a process want to remove this
	static void RemoveProcess(std::string ID)
	{
		std::map<std::string, timespec>::iterator Location = pData.begin();
		Location = pData.find(ID);
		if(Location != pData.end())
		{
			pData.erase(Location);
		}
	}

	static void StartProcessTimer(std::string ID)
	{
		std::map<std::string, timespec>::iterator Location = StartTime.begin();
		Location = StartTime.find(ID);
		if(Location != StartTime.end())
		{
			timespec Start_Time;
			clock_gettime(CLOCK_REALTIME, &Start_Time);
			Location->second.tv_nsec = Start_Time.tv_nsec;
			Location->second.tv_sec = Start_Time.tv_sec;
		}
	}

	static void EndProcessTimer(std::string ID)
	{
		timespec Time_End, Time_Total;
		clock_gettime(CLOCK_REALTIME, &Time_End);

		std::map<std::string, timespec>::iterator Location = StartTime.begin();
		Location = StartTime.find(ID);
		if(Location != StartTime.end())
		{
			Time_Total.tv_sec = Time_End.tv_sec - Location->second.tv_sec;
			Time_Total.tv_nsec = Time_End.tv_nsec - Location->second.tv_nsec;
			SetProcessDuration(Time_Total, ID);
		}
	}

	static void SetProcessDuration(timespec t, std::string ID)
	{
		//Find the process entry we are looking for
		std::map<std::string, timespec>::iterator Location = pData.begin();
		Location = pData.find(ID);
		if(Location != pData.end())
		{
			Location->second.tv_nsec += t.tv_nsec;
			Location->second.tv_sec += t.tv_sec;
			if(Location->second.tv_nsec >= 1e9)
			{
				Location->second.tv_sec += 1;
				Location->second.tv_nsec -= 1e9;
			}
		}
	}

	static void ClearTime()
	{
		std::map<std::string, timespec>::iterator Location = pData.begin();
		if(Location != pData.end())
		{
			Location->second.tv_nsec = 0;
			Location->second.tv_sec = 0;
		}
	}

	static void GetProfileData()
	{
		if(!IsEnabled())
		{
			std::cout << "Profiling disabled in libmerlin. Build Merlin with -DMERLIN_PROFILE to enable" << std::endl;
		}
		double total_t = 0.0;
		std::cout << std::endl << std::setw(24) << std::left << "PROCESS NAME\tTIME" << std::endl;
		std::map<std::string, timespec>::iterator Location = pData.begin();
		Location = pData.find("TOTAL");
		if(Location != pData.end())
		{
			total_t = ((double) Location->second.tv_nsec / 1.0e9) + (double) Location->second.tv_sec;
		}
		Location = pData.begin();
		while(Location != pData.end())
		{
			double time = ((double) Location->second.tv_nsec / 1.0e9) + (double) Location->second.tv_sec;
			std::cout << std::setw(24) << std::left << Location->first << "\t" << time << "\t" << time * 100.0
				/ total_t << "%" << std::endl;
			Location++;
		}
	}

	static bool IsEnabled();

private:
	static std::map<std::string, timespec> pData;
	static std::map<std::string, timespec> StartTime;

};

#else
class MerlinProfile
{
public:
	static void AddProcess(std::string ID)
	{
	}
	static void RemoveProcess(std::string ID)
	{
	}
	static void StartProcessTimer(std::string ID)
	{
	}
	static void EndProcessTimer(std::string ID)
	{
	}
	static void SetProcessDuration(timespec t, std::string ID)
	{
	}
	static void ClearTime()
	{
	}
	static void GetProfileData()
	{
		std::cout << "Profiling disabled. Build Merlin with -DMERLIN_PROFILE to enable" << std::endl;
	}
	static bool IsEnabled();
private:
	static std::map<std::string, timespec> pData;
	static std::map<std::string, timespec> StartTime;
};

#endif
#endif
