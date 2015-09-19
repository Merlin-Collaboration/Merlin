#ifndef _APERTURECONFIGURATION_H_
#define _APERTURECONFIGURATION_H_

#include <vector>
#include <string>
#include <fstream>

#include "AcceleratorModel/AcceleratorModel.h"

using namespace std;

class ApertureConfiguration
{
public:

	struct ap
	{
		double s;
		double ap1;
		double ap2;
		double ap3;
		double ap4;
	};

	ap ApertureEntry;
	std::vector<ap> ApertureList;
	std::vector<ap> ApertureListForElement;

	ApertureConfiguration();
	ApertureConfiguration(string);
	ApertureConfiguration(string, bool);
	void LoadApertureConfiguration(string);
	void ConfigureElementApertures(AcceleratorModel*);
	
	void SetAllRectEllipse(){allRectEllipse=1;}

	string ApertureType;
	//Output log
	ostream* log;
	bool logFlag;
	
	bool allRectEllipse;

	//Set the stream for the collimator settings log.
	void SetLogFile (ostream& os);

	//Enable/disable logging
	void EnableLogging(bool);
};

#endif
