//----------------------------------------------
//    How to use coupler wakefields
//    DK 29.11.2008
//

#include "Random/RandomNG.h"
#include "BeamModel/BeamData.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/Components.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"

#include "TeslaCoupler.h"
#include "BeamDynamics/ParticleTracking/CouplerWakeFieldProcess.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"

#include <fstream>
#include <iostream>

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;


AcceleratorModel* ConstructLinac(WakePotentials*,BeamData&);

int main() {

  RandomNG::init();
  
//-----------------------------------------------------
//             define the initial beam
//             at 15 GeV
//-----------------------------------------------------

  double E0    = 15*GeV;
  double gamma = E0/MeV/ElectronMassMeV;

  BeamData beam;
  beam.p0      = E0;
  beam.charge  = 2.0e10;
  beam.emit_x  = 8.0e-06/gamma;
  beam.emit_y  = 0.02e-06/gamma;
  beam.sig_dp  = 0.0107;
  beam.sig_z   = 300.0e-06;

  int npart = 10000;

//----------------------------------------------
//             construct the accelerator model,
//             a simple ILC like linac 
//----------------------------------------------

  TeslaCoupler*   wake  = new TeslaCoupler;
  AcceleratorModel* model = ConstructLinac(wake,beam);

  CouplerWakeFieldProcess* couplerWakeProc = new CouplerWakeFieldProcess(1);
  couplerWakeProc->ApplyImpulseAt(ParticleTracking::CouplerWakeFieldProcess::atExit);

  // set responsibility
  wake->SetExpectedProcess(couplerWakeProc);
  
//-----------------------------------------------------
//             Construct bunch
//-----------------------------------------------------

  ParticleBunchConstructor* PBC = new ParticleBunchConstructor(beam,npart,normalDistribution);
  PBC->ForceCentroid(true); // set centroid to 0
  ParticleBunch* startBunch1 = PBC->ConstructParticleBunch();
  ParticleBunch* startBunch2 = PBC->ConstructParticleBunch();

  PSmoments S;
  startBunch1->GetMoments(S);
        
  double E = startBunch1->GetReferenceMomentum();
         E*=1+S.mean(ps_DP);
  double gam = E/MeV/ElectronMassMeV;
  double emitx=sqrt(S(0,0)*S(1,1)-pow(S(0,1),2));
  double emity=sqrt(S(2,2)*S(3,3)-pow(S(2,3),2));
  cout<<"initial energy "<<E<<endl;
  cout<<"initial gamma*emittance_x "<<emitx*gam<<endl;
  cout<<"initial gamma*emittance_y "<<emity*gam<<endl;

//----------------------------------------------
//             Construct Tracker and add
//             the combined coupler and cavity 
//             wakefield processes
//----------------------------------------------

  ParticleTracker* tracker = new ParticleTracker(model->GetBeamline());

//----------------------------------------------
//  Final bunch parameters without wakefields
//----------------------------------------------
   
  ParticleBunch* finalBunch1 = tracker->Track(startBunch1);

  PSmoments S2;
  finalBunch1->GetMoments(S2);
        
  E = finalBunch1->GetReferenceMomentum();
  E*=1+S2.mean(ps_DP);
  gam = E/MeV/ElectronMassMeV;
  emitx=sqrt(S2(0,0)*S2(1,1)-pow(S2(0,1),2));
  emity=sqrt(S2(2,2)*S2(3,3)-pow(S2(2,3),2));
  cout<<"No wakefields "<<endl;
  cout<<"final energy "<<E<<endl;
  cout<<"final gamma*emittance_x "<<emitx*gam<<endl;
  cout<<"final gamma*emittance_y "<<emity*gam<<endl;


//----------------------------------------------
// Final bunch parameters with wakefields
//----------------------------------------------

  //now we track with coupler wakefields included
  tracker->AddProcess(couplerWakeProc);
  ParticleBunch* finalBunch2 = tracker->Track(startBunch2);
  
  PSmoments S3;
  finalBunch2->GetMoments(S3);
  E = finalBunch2->GetReferenceMomentum();
  E*=1+S3.mean(ps_DP);
  gam = E/MeV/ElectronMassMeV;
  emitx=sqrt(S3(0,0)*S3(1,1)-pow(S3(0,1),2));
  emity=sqrt(S3(2,2)*S3(3,3)-pow(S3(2,3),2));
  cout<<"With cavity and (overestimated) strong coupler wakefields "<<endl;
  cout<<"final energy "<<E<<endl;
  cout<<"final gamma*emittance_x "<<emitx*gam<<endl;
  cout<<"final gamma*emittance_y "<<emity*gam<<endl;
  
  return 0;

}
//----------------------------------------------
// Create a simple ILC like linac
//   300 quads - 33.9 m distance 
//   24 cavities between each quad pair
AcceleratorModel* ConstructLinac(WakePotentials* wake, BeamData& beam){

  AcceleratorModelConstructor* accel_modConst = new AcceleratorModelConstructor();
  accel_modConst->NewModel();

  // drift between cavities
  double driftLen = 0.3346*meter;

  // cavity parameters
  double cavLen =   1.0362*meter;
  double volt   =  32.6*MV;
  double phase  =  -0.0925025;
  double freq   =   1.3e+09;
  double eloss  =   0.0458*Volt;

  // quadrupole parameters
  double quadLen =  0.666*meter;
  double f_k1    =  0.0524;
  double d_k1    = -0.0471;
  double E       = beam.p0;

  double L=0;

  for(int q=0;q<302;q++){
  
     double brho=E/eV/SpeedOfLight;

     if(q%2==0){
        accel_modConst->AppendComponent(new Quadrupole("FQ",quadLen,f_k1*brho) );
        //cout<<"FQ "<<L<<" "<<E<<endl;
     } else {
        accel_modConst->AppendComponent(new Quadrupole("DQ",quadLen,d_k1*brho) );
        //cout<<"DQ "<<L<<" "<<E<<endl;
     }
     L+=quadLen;
     
     for(int c=0;c<24;c++){
        accel_modConst->AppendComponent(new Drift("DriftCav", driftLen) );
        TWRFStructure* cavity = 
          accel_modConst->AppendComponent(new TWRFStructure("Cav",cavLen,freq,volt/cavLen,phase) );
        cavity->SetWakePotentials(wake);
        E+=cos(phase)*volt-eloss;
        L+=driftLen+cavLen;
     }
     accel_modConst->AppendComponent(new Drift("DriftCav", driftLen) );
     L+=driftLen;

  }
  
  // set twiss parameters for this lattice
  beam.beta_x  = 106.37*meter;
  beam.beta_y  =  37.27*meter;
  beam.alpha_x = -1.857;
  beam.alpha_y =  0.6584;

  return accel_modConst->GetModel();

}
