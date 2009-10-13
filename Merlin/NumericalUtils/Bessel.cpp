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
#include <stdlib.h>

// Error function and complimentry error function
// Taken for NRiC (and modified)

double BesselI0(double x)
{
    double besselI0 = 0;;
    static const double p1 = 1.0;
    static const double p2 = 3.5156229;
    static const double p3 = 3.0899424;
    static const double p4 = 1.2067492;
    static const double p5 = 0.2659732;
    static const double p6 = 0.0360768;
    static const double p7 = 0.0045813;
    static const double q1 = 0.39894228;
    static const double q2 = 0.01328592;
    static const double q3 = 0.00225319;
    static const double q4 =-0.00157565;
    static const double q5 = 0.00916281;
    static const double q6 =-0.02057706;
    static const double q7 = 0.02635537;
    static const double q8 =-0.01647633;
    static const double q9 = 0.00392377;

    double ax = fabs(x);
    if(ax<3.75)
    {
        double y = x*x/14.0625;
        besselI0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
    }
    else
    {
        double y = 3.75/ax;
        besselI0 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
    }
    return besselI0;
}

double BesselI1(double x)
{
    double besselI1 = 0;;
    static const double p1 = 0.5;
    static const double p2 = 0.87890594;
    static const double p3 = 0.51498869;
    static const double p4 = 0.15084934;
    static const double p5 = 0.02658733;
    static const double p6 = 0.00301532;
    static const double p7 = 0.00032411;
    static const double q1 = 0.39894228;
    static const double q2 =-0.03988024;
    static const double q3 =-0.00362018;
    static const double q4 = 0.00163801;
    static const double q5 =-0.01031555;
    static const double q6 = 0.02282967;
    static const double q7 =-0.02895312;
    static const double q8 = 0.01787654;
    static const double q9 =-0.00420059;

    double ax = fabs(x);
    if(ax<3.75)
    {
        double y = x*x/14.0625;
        besselI1 = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
    }
    else
    {
        double y = 3.75/ax;
        besselI1 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
        if(x<0)
            besselI1 *= -1.;
    }
    return besselI1;
}

double BesselIn(int n, double x)
{
#define IACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

    if(n==0)
        return BesselI0(x);
    if(x==0)
        return 0;
    if(n==1)
        return BesselI1(x);

    double besselI = 0;
    double tox = 2.0/fabs(x);
    double bip = 0.0;
    double bi  = 1.0;
    int m = 2*(n+(int)sqrt(IACC*n));
    for(int j=m; j>0; j--)
    {
        double bim = bip + j*tox*bi;
        bip = bi;
        bi  = bim;
        if(fabs(bi)>BIGNO)
        {
            besselI *= BIGNI;
            bi  *= BIGNI;
            bip *= BIGNI;
        }
        if(j==n)
            besselI=bip;
    }
    besselI *= BesselI0(x)/bi;
    if(x<0 && div(n,2).rem==1)
        besselI *= -1.;

    return besselI;
}
