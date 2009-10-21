
#include "BeamModel/BeamData.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "SpoilerWakeProcess.h"
//#include "BeamDynamics/ParticleTracking/SpoilerWakeProcess.h"
#include "Random/RandomNG.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/WakePotentials.h"
#include "AcceleratorModel/SpoilerWakePotentials.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/StdComponent/Drift.h"
//#include "TaperedSpoilerWakePotentials.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"
#include <typeinfo>
#include <iostream>
using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

//------------------------------------------------------------------------------------------------------------
//      the resistive wake potentials  (in MKS ssytem)
//------------------------------------------------------------------------------------------------------------

class ResistivePotential:public SpoilerWakePotentials  {
public:

    double sigma,b,leng,scale,step;
    int ncoeff;
    double lo,hi,* coeff;
    ResistivePotential(int m, double ss, double bb, double l,char*filename):SpoilerWakePotentials(m, 0., 0.){
       sigma=ss;    b=bb;
   double Z0=377;
   // double fourpieps=4*3.1416*8.854E-12;
   double fourpieps=1; // pu in later
   scale=pow(2*b*b/(Z0*sigma),1./3.);
   leng=l;
   cout<<" resistive collimator radius "<<b<<" length "<<leng<<" conductivity "<<sigma<<endl;
   cout<<" Scale length "<<scale<<endl;
   ifstream file;
   file.open(filename);
   file>>ncoeff;
   file>>lo;
   file>>hi;
   step=(hi-lo)/(ncoeff-1);
   cout<<" Resistive wall table has "<<ncoeff<<" values from "<<lo<<" to "<<hi<<endl;
   coeff=new double[ncoeff];
   for(int i=0;i<ncoeff;i++) file>>coeff[i];
      }       
    double Wlong (double z) const {return 1;};
    double Wtrans (double z) const { return 1;};
    double Wtrans (double z,int m) const {if(z<0) return 0;
   double c=3.E8;
   double Pi=3.14159;
   double Z0=377;
   double s=z/scale;
   double E;
   //double fourpieps=4*3.1416*8.854E-12;
   double fourpieps=1; // pu in later
      //double Chao=scale*(1/(Pi*b))*sqrt(2*Pi*scale/(b*b))*leng/sqrt(z);
      double Chao=2*leng*scale*(1/(fourpieps*sqrt(2*pi)))*sqrt(scale/z)/(b*b);
     int index=int(0.5+(s-lo)/step);
     if(index<0) index=0;
     int i1;
     if(index==0){i1=1;} else if(index==ncoeff-1) {i1=ncoeff-2;} else {i1=index;}
     if(index>=ncoeff) {E=Chao;} else 
         {
          double d=(s-(lo+i1*step))/step;
          E=coeff[i1]+d*(coeff[i1+1]-coeff[i1-1])/2+
                 d*d*(coeff[i1+1]+coeff[i1-1]-2*coeff[i1])/2;
           E=-scale*leng*E; // minus sign should have been in Mathematica
           E=E/fourpieps; // minus sign should have been in Mathematica
           E=E/(b*b);
             }
     cout<<" kick "<<z<<" "<<E<<" "<<Chao<<endl;
     E=2*E/(b*b); // naive way to get from mode 0 to mode 1
     return E;};
    double Wlong (double z,int m) const { return z>0? 1:0;};
};



class ResistiveWakePotentials:public SpoilerWakePotentials {
public: 
  
   double* coeff;
   double rad, sigma, length;
  
   ResistiveWakePotentials(int m, double r, double s, double l):SpoilerWakePotentials(m, r, s)
   {
    rad = r;
    sigma = s;
    length = l;
    cout<<" Making new ResistiveWakePotentials with length "<<length<<endl;
    coeff = new double[m+1];
 
    int delta;
     if (m==0)
         delta=1;
     else
         delta = 0;

     for (int i=0; i<(m+1); i++)
        {
         coeff[i]= 1/(pi*pow(rad,2*i+1)*(1+delta));
         }
   }

   double Wlong (double z) const {return 1;};
   double Wtrans (double z) const { return 1;};
   double Wtrans (double z, int m) const {cout<<" call for "<<m<<endl; 
    return z>0? coeff[m]*sqrt(376.74/(pi*sigma))*length/sqrt(z):0;};
   double Wlong (double z, int m) const {return z>0? -coeff[m]*sqrt(1/sigma*376.6)*sqrt(z)*length:0;};

};


