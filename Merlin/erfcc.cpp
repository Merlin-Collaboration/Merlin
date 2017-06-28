/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
//
/////////////////////////////////////////////////////////////////////////

#include "utils.h"
#include <cmath>

using namespace std;

double NormalBin(double x1, double x2)
{
	static const double root2 = sqrt(2.0);
	return 0.5*(erfc(x1/root2)-erfc(x2/root2));
}
