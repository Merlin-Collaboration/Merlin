//--------------
// AcceleratorErrors
// DK 1.12.08

#ifndef AcceleratorErrors_h
#define AcceleratorErrors_h 1

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Frames/LatticeFrame.h"
#include "AcceleratorModel/Frames/ComponentFrame.h"
#include "AcceleratorModel/Frames/SequenceFrame.h"
#include "AcceleratorModel/Supports/SupportStructure.h"

#include "utility/StringPattern.h"
#include "NumericalUtils/NumericalConstants.h"
#include "Random/RandomNG.h"

#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

class AcceleratorErrors {
public:
	AcceleratorErrors(){};
	void SetErrors(double  xrms=0, double  yrms=0, double  zrms=0,
	          double meanx=0, double meany=0, double meanz=0 ){
		vx=xrms*xrms; vy=yrms*yrms; vz=zrms*zrms;
		mx=meanx; my=meany; mz=meanz;
		clear=true;
	}
	void AddErrors(double  xrms=0, double  yrms=0, double  zrms=0,
	          double meanx=0, double meany=0, double meanz=0 ){
		vx=xrms*xrms; vy=yrms*yrms; vz=zrms*zrms;
		mx=meanx; my=meany; mz=meanz;
		clear=false;
	}
	void ApplyShifts(AcceleratorModel::Beamline& b, const string& p);
	void ApplyRotations(AcceleratorModel::Beamline& b, const string& p);
	void Report(ostream* l){log=l;};
private:
        bool      clear;
        double vx,vy,vz;
        double mx,my,mz;
	ostream* log;
};
#endif
