#ifndef _MerlinProfile_hpp_
#define _MerlinProfile_hpp_ 1

#include <ctime>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
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
		pData.insert(std::pair<std::string,timespec>(ID,pTime));
		StartTime.insert(std::pair<std::string,timespec>(ID,pTime));
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
	//AddProcess(ID);
	std::map<std::string, timespec>::iterator Location = StartTime.begin();
	Location = StartTime.find(ID);
	if(Location != StartTime.end())
	{
		//std::cout << "Adding Profile for: " << ID << std::endl;
		timespec Start_Time;
		clock_gettime(CLOCK_REALTIME, &Start_Time);
		//StartTime.insert(std::pair<std::string,timespec>(ID,Start_Time));
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
	        Time_Total.tv_sec = Time_End.tv_sec - Location->second.tv_sec;;
	        Time_Total.tv_nsec = Time_End.tv_nsec - Location->second.tv_nsec;;
		SetProcessDuration(Time_Total,ID);
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
	double total_t = 0.0;
	std::cout << std::endl << std::setw(24) << std::left << "PROCESS NAME\tTIME" << std::endl;
	std::map<std::string, timespec>::iterator Location = pData.begin();
	Location = pData.find("TOTAL");
	if(Location != pData.end())
	{
		total_t = ((double)Location->second.tv_nsec/1.0e9) + (double)Location->second.tv_sec;
	}
	Location = pData.begin();
	while(Location !=pData.end())
	{
		double time = ((double)Location->second.tv_nsec/1.0e9) + (double)Location->second.tv_sec;
		std::cout << std::setw(24) << std::left << Location->first << "\t" << time << "\t" << time*100.0/total_t << "%" << std::endl;		
		Location++;
	}
}
private:
static std::map<std::string, timespec> pData;
static std::map<std::string, timespec> StartTime;

};
//std::map<std::string, timespec> MerlinProfile::pData;
//std::map<std::string, timespec> MerlinProfile::StartTime;
#endif
