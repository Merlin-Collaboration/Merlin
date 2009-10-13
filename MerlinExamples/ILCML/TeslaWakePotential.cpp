#include <cmath>

#include "TeslaWakePotential.h"

//using namespace std;



double TeslaWakePotentials::Wlong(double z) const

{

	return 38.1e+12*(1.165*exp(-sqrt(z/3.65e-3))-0.165);

}



double TeslaWakePotentials::Wtrans(double z) const
{
	// Original TDR transverse wake
	return 1.0e+12*(1290.0*sqrt(z)-2600.0*z);
	// New transverse wake
	//	double arZ=sqrt(z/0.92e-03);
	//	return 1.21e14*(1.0-(1.0+arZ)*exp(-arZ));
}

