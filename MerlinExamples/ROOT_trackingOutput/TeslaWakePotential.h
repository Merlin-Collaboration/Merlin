// class TeslaWakePotential

//

// A form of WakeFieldProcess::WakePotential that uses linear interpolation

// from a tabulated wakefield file.

//



#ifndef _H_TeslaWakePotential

#define _H_TeslaWakePotential



#include "AcceleratorModel/WakePotentials.h"



class TeslaWakePotentials : public WakePotentials {

public:

	double Wlong(double z) const;

	double Wtrans(double z) const;

};





#endif



