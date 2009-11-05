
#include "BeamModel/BeamData.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "Collimators/SpoilerWakeProcess.h"
#include "Random/RandomNG.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/WakePotentials.h"
#include "Collimators/SpoilerWakePotentials.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/StdComponent/Drift.h"
//#include "TaperedSpoilerWakePotentials.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"
#include "collimatortable.h"
#include <typeinfo>
#include <iostream>
#include <strstream>
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
    collimatortable** Transverse;
    collimatortable** Longitudinal;
    ResistivePotential(int m, double ss, double bb, double l,char*filename,double tau=0):SpoilerWakePotentials(m, 0., 0.){
       sigma=ss;    b=bb;
       double Z0=377;
       scale=pow(2*b*b/(Z0*sigma),1./3.);
       leng=l;
       cout<<" resistive collimator radius "<<b<<" length "<<leng<<" conductivity "<<sigma<<" Scale length "<<scale;
       double xi=pow(scale/b,2);
       double Gamma=SpeedOfLight*tau/scale;
       cout<<" xi "<<xi<<" Gamma "<<Gamma<<endl; 
       Transverse=new collimatortable*[m+1];
       Longitudinal=new collimatortable*[m+1];
          for(int mode=0;mode<=m;mode++){
	   char line[100]={100*0};
           ostrstream ss(line,100);
           ss<<filename<<"L"<<mode<<".txt";
           Longitudinal[mode]=new collimatortable(line,Gamma,xi);
           if(mode>0) {
	      char line[100]={100*0};
              ostrstream ss(line,100);
              ss<<filename<<"T2m"<<mode<<".txt";
              Transverse[mode]=new collimatortable(line,Gamma,xi);
              }
	      }
       }// end of constructor
    double Wlong (double z) const {return 1;};
    double Wtrans (double z) const { return 1;};
    double Wtrans (double z,int m) const {
     if(z<0) return 0;
     double s=z/scale;
     double fourpieps=1; // pu in later
     double Chao=-2*(1/(sqrt(2*pi)))*sqrt(scale/z);
     double E=Transverse[m]->inrange(s) ? Transverse[m]->interpolate(s) : Chao; 
    // cout<<m<<" "<<s<<" "<<E1<<" "<<Chao<<endl;
     E= scale*leng*E; // minus sign should have been in Mathematica
     E=E/fourpieps; // minus sign should have been in Mathematica
     E=E/pow(b,2*m+2); 
    // cout<<m<<" "<<s<<" "<<E<<endl;
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


