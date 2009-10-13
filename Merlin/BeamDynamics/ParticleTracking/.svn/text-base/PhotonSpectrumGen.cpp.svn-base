// GeSynRad.C, GEnerator for SYNchrotron RADiation, distribution version in C++
// by Helmut Burkhardt, CERN
// there is also a FORTRAN equivalent (dsynrad.f)
//
// reference:
// H.Burkhardt, Monte Carlo generator for synchrotron radiation
// CERN-LEP-Note 632 (Dec 1990)
// available from http://wwwcn.cern.ch/~hbu/Welcome.html
//
// Modified by N. Walker (DESY) to use with the Merlin Class Library
//

#include <cmath>
#include "NumericalUtils/NumericalConstants.h"
#include "Random/ACG.h"

namespace {
// HB spectrum, Chebyshev series
double SynRadC(double x);
double SynGenC(double x);

// Random number generator [0-1]
double Ran1();
};

double HBSpectrumGen(double uc)
{
    return uc*SynGenC(0.0);
}

double AWSpectrumGen(double uc)
{
    double u1 = Ran1();
    double u2 = 0.57 * uc * pow(u1,-1./3) * pow((1-u1),pi);
    return u2;
}

namespace {

double Ran1()
{
    static ACG pgen(12345);
    return pgen.asDouble();
}


double SynGenC(double xmin)
// C++ version
// synchrotron radiation photon spectrum generator
// returns photon energy in units of the critical energy
// xmin is the lower limit for the photon energy to be generated
//          see H.Burkhardt, LEP Note 632
{
    static bool DoInit=true;
    static double a1,a2,c1,xlow,ratio,LastXmin=-1.;

    if(DoInit || xmin!=LastXmin) {
        DoInit=false;
        LastXmin=xmin;
        if(xmin>1.)
            xlow=xmin;
        else
            xlow=1.;

        // initialize constants used in the approximate expressions
        // for SYNRAD   (integral over the modified Bessel function K5/3)
        a1=SynRadC(1.e-38)/pow(1.e-38,-2./3.); // = 2**2/3 GAMMA(2/3)
        a2=SynRadC(xlow)/exp(-xlow);
        c1=pow(xmin,1./3.);
        // calculate the integrals of the approximate expressions
        if(xmin<1.) { // low and high approx needed
            double sum1ap=3.*a1*(1.-pow(xmin,1./3.)); //     integral xmin --> 1
            double sum2ap=a2*exp(-1.);                    //     integral 1 --> infin
            ratio=sum1ap/(sum1ap+sum2ap);
        }
        else {// only high approx needed
            ratio=0.;     //     generate only high energies using approx. 2

            a2=SynRadC(xlow)/exp(-xlow);
        }
    }
    // Init done, now generate
    double appr,exact,result;
    do {
        if(Ran1()<ratio) {          // use low energy approximation
            result=c1+(1.-c1)*Ran1();
            result*=result*result;  // take to 3rd power;
            exact=SynRadC(result);
            appr=a1*pow(result,-2./3.);
        }
        else {                     // use high energy approximation
            result=xlow-log(Ran1());
            exact=SynRadC(result);
            appr=a2*exp(-result);
        }
    }
    while(exact<appr*Ran1());    // reject in proportion of approx
    return result;               // result now exact spectrum with unity weight
}

//	x :    energy normalized to the critical energy
//	returns function value SynRadC   photon spectrum dn/dx
//	(integral of modified 1/3 order Bessel function)
//	principal: Chebyshev series see H.H.Umstaetter CERN/PS/SM/81-13 10-3-1981
//	see also my LEP Note 632 of 12/1990
//	converted to C++, H.Burkhardt 21-4-1996    */

double SynRadC(double x)
{
    double synrad=0.;
    if((x>0.)&&(x<800.)) { // otherwise result synrad remains 0
        if(x<6.) {
            double a,b,z;
            z=x*x/16.-2.;
            b=          .00000000000000000012;
            a=z*b  +    .00000000000000000460;
            b=z*a-b+    .00000000000000031738;
            a=z*b-a+    .00000000000002004426;
            b=z*a-b+    .00000000000111455474;
            a=z*b-a+    .00000000005407460944;
            b=z*a-b+    .00000000226722011790;
            a=z*b-a+    .00000008125130371644;
            b=z*a-b+    .00000245751373955212;
            a=z*b-a+    .00006181256113829740;
            b=z*a-b+    .00127066381953661690;
            a=z*b-a+    .02091216799114667278;
            b=z*a-b+    .26880346058164526514;
            a=z*b-a+   2.61902183794862213818;
            b=z*a-b+  18.65250896865416256398;
            a=z*b-a+  92.95232665922707542088;
            b=z*a-b+ 308.15919413131586030542;
            a=z*b-a+ 644.86979658236221700714;
            double p;
            p=.5*z*a-b+  414.56543648832546975110;
            a=          .00000000000000000004;
            b=z*a+      .00000000000000000289;
            a=z*b-a+    .00000000000000019786;
            b=z*a-b+    .00000000000001196168;
            a=z*b-a+    .00000000000063427729;
            b=z*a-b+    .00000000002923635681;
            a=z*b-a+    .00000000115951672806;
            b=z*a-b+    .00000003910314748244;
            a=z*b-a+    .00000110599584794379;
            b=z*a-b+    .00002581451439721298;
            a=z*b-a+    .00048768692916240683;
            b=z*a-b+    .00728456195503504923;
            a=z*b-a+    .08357935463720537773;
            b=z*a-b+    .71031361199218887514;
            a=z*b-a+   4.26780261265492264837;
            b=z*a-b+  17.05540785795221885751;
            a=z*b-a+  41.83903486779678800040;
            double q;
            q=.5*z*a-b+28.41787374362784178164;
            double y;
            double const twothird=2./3.;
            y=pow(x,twothird);
            synrad=(p/y-q*y-1.)*1.81379936423421784215530788143;
        }
        else {// 6 < x < 174
            double a,b,z;
            z=20./x-2.;
            a=      .00000000000000000001;
            b=z*a  -.00000000000000000002;
            a=z*b-a+.00000000000000000006;
            b=z*a-b-.00000000000000000020;
            a=z*b-a+.00000000000000000066;
            b=z*a-b-.00000000000000000216;
            a=z*b-a+.00000000000000000721;
            b=z*a-b-.00000000000000002443;
            a=z*b-a+.00000000000000008441;
            b=z*a-b-.00000000000000029752;
            a=z*b-a+.00000000000000107116;
            b=z*a-b-.00000000000000394564;
            a=z*b-a+.00000000000001489474;
            b=z*a-b-.00000000000005773537;
            a=z*b-a+.00000000000023030657;
            b=z*a-b-.00000000000094784973;
            a=z*b-a+.00000000000403683207;
            b=z*a-b-.00000000001785432348;
            a=z*b-a+.00000000008235329314;
            b=z*a-b-.00000000039817923621;
            a=z*b-a+.00000000203088939238;
            b=z*a-b-.00000001101482369622;
            a=z*b-a+.00000006418902302372;
            b=z*a-b-.00000040756144386809;
            a=z*b-a+.00000287536465397527;
            b=z*a-b-.00002321251614543524;
            a=z*b-a+.00022505317277986004;
            b=z*a-b-.00287636803664026799;
            a=z*b-a+.06239591359332750793;
            double p;
            p=.5*z*a-b    +1.06552390798340693166;
            double const pihalf=pi/2.;
            synrad=p*sqrt(pihalf/x)/exp(x);
        }
    }
    return synrad;
}

};

