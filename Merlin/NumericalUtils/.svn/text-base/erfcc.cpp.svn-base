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

#include "NumericalUtils/utils.h"
#include <cmath>

using namespace std;

// Error function and complimentry error function
// Taken for NRiC (and modified)

double erfc(double x)
{
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                                 t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                                                                 t*(-0.82215223+t*0.17087277)))))))));
    return  x >= 0.0 ? ans : 2.0-ans;
}

double NormalBin(double x1, double x2)
{
    static const double root2 = sqrt(2.0);
    return 0.5*(erfc(x1/root2)-erfc(x2/root2));
}
