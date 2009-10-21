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

class TaperedCollimatorPotentials:public SpoilerWakePotentials  {
public:
    double a,b;
    double* coeff;
    TaperedCollimatorPotentials(int m, double rada, double radb):SpoilerWakePotentials(m,0.,0.){
      a=rada;
      b=radb;
      coeff=new double[m+1];
      for(int i=0;i<=m;i++) {coeff[i]=2*(1./pow(a,2*i)-1./pow(b,2*i));}      
     }
    ~TaperedCollimatorPotentials(){ delete[] coeff;}
    double Wlong (double z) const {return 0;};
    double Wtrans (double z) const {return 0;};
    double Wlong (double z,int m) const {return  z>0?-(m/a)*coeff[m]/exp(m*z/a):0;};
    double Wtrans (double z,int m) const { return z>0?coeff[m]/exp(m*z/a):0;};
};

