//----------------------------------------------
//    example on using collimator wakefields together with cavity wakefield 
//    DK 29.11.2008
//
#include "BeamModel/BeamData.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/SpoilerWakeProcess.h" 
#include "BeamDynamics/ParticleTracking/CollimateParticleProcess.h"
#include "Random/RandomNG.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/SpoilerWakePotentials.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"

#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "TeslaWakePotential.h"

#include <fstream>
#include <iostream>

// these classes implement examples of wake potentials
// see A.M. Toader et al., EPAC08, Genua, WEPP161 for 
// a collection of collimator wakefield formulae 
#include "AcceleratorModel/SpoilerPotentialModels.h"
  
using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

int main() {

  RandomNG::init();
  
  // initial offset in meters
  double offset = 0.0015;
  // number of modes
  int modes = 5;
   
//-----------------------------------------------------
//             define the beam
//             SLAC test with 1.19 GeV
//-----------------------------------------------------

  BeamData mybeam;
  mybeam.charge = 2.0e10;
  mybeam.beta_x = 3.*meter;  
  mybeam.beta_y = 10.*meter;  
  mybeam.emit_x = 0.36*millimeter; 
  mybeam.emit_y = 0.16*millimeter; 
  mybeam.sig_z  = 0.65*millimeter;        
  mybeam.p0     = 1.19*GeV;
     
  mybeam.y0 = offset;

  int npart = 10000;

//-----------------------------------------------------
//             Construct bunch
//-----------------------------------------------------          
  ParticleBunchConstructor* PBC = new ParticleBunchConstructor(mybeam,npart,normalDistribution);
  PBC->ForceCentroid(true); // set centroid to 0
  ParticleBunch* startBunch = PBC->ConstructParticleBunch();

//-----------------------------------------------------
// the initial mean and sigma values
//-----------------------------------------------------

  double xmeani = startBunch->GetMoments(0).first;
  double xsigi  = startBunch->GetMoments(0).second;
  double ymeani = startBunch->GetMoments(2).first;
  double ysigi  = startBunch->GetMoments(2).second;

  cout<<"Initial bunch parameters:"<<endl;
  cout<<"(beta_x= "<<mybeam.beta_x<<", emit_x="<<mybeam.emit_x<<" => sig_x="
      <<sqrt(mybeam.beta_x*mybeam.emit_x)<<" [m]"<<endl;
  cout<<" beta_y="<<mybeam.beta_y<<", emit_y="<<mybeam.emit_y<<" => sig_y="
      <<sqrt(mybeam.beta_y*mybeam.emit_y)<<" [m]"<<endl;
  cout<<" offset= "<<offset<<" [m])"<<endl<<endl;
  cout<<"Mean x ="<<xmeani<<' '<<"Sigma x ="<<xsigi<<endl;
  cout<<"Mean y ="<<ymeani<<' '<<"Sigma y ="<<ysigi<<endl;
  cout<<"yp angle :"<<startBunch->GetMoments(3).first<<endl<<endl;;

//-----------------------------------------------------
//             construct the accelerator model
//             a spoiler between 2 drifts
//             + 1 cavity with Tesla wake potential
//-----------------------------------------------------

  AcceleratorModelConstructor* accelerator_model = new AcceleratorModelConstructor();
  accelerator_model->NewModel();

  double driftlength1 = 1.0*meter;
  Drift* drift1       = new Drift("aDrift1", driftlength1);

  double spoilerlength  = 177.*millimeter;  
  double spoilerthick   = 0.004*meter;  
  Spoiler* spoiler      = new Spoiler("aSpoiler", spoilerlength, spoilerthick);
  double aperturewidth  = 1.9*millimeter;  
  double apertureheight = 1.9*millimeter; 
  RectangularAperture* aperture = new RectangularAperture(aperturewidth, apertureheight);
  spoiler->SetAperture(aperture);

  double driftlength2 = 1.0*meter;
  Drift* drift2       = new Drift("aDrift2", driftlength2);

  // we add just one Tesla cavity as an example
  double len   = 1.0362;
  double volt  = 0.0315665;
  double phase = -0.0925025;
  double freq  = 1.3e+09;
  TWRFStructure* cavity = new TWRFStructure("mycav",len,freq,volt/len,phase);

  accelerator_model->AppendComponent(cavity);
  accelerator_model->AppendComponent(drift1);
  accelerator_model->AppendComponent(spoiler);
  accelerator_model->AppendComponent(drift2);

  // the spoiler needs a CollimateParticleProcess
  // and we want to have scattering 
  ofstream lossummary("loss_summary.dat");
  CollimateParticleProcess* collimation = new CollimateParticleProcess(0,COLL_AT_EXIT, &lossummary);
  collimation->ScatterAtSpoiler(true); 

  AcceleratorModel* model = accelerator_model->GetModel();
 
//-----------------------------------------------------
//             Create a wake potential and
//             add it to the spoiler element
//-----------------------------------------------------
      
 // apply the geometric wakefields 
 TaperedCollimatorPotentials* collWake 
     =  new TaperedCollimatorPotentials(modes,aperturewidth/2,apertureheight/2);
 spoiler->SetWakePotentials(collWake);
    
 // //another example: a resistive wakefield
 // double conductivity = 2.38e6;         //titanium 
 // double conductivity = 5.98*10000000;  //copper
 // double conductivity = 3.08e7;         //berrilium
 // double conductivity = 6.e4;           //carbon 
 // double conductivity = 4.5e6;          //for TiN
 // ResistiveWakePotentials* resWake =  new ResistiveWakePotentials(modes, aperturewidth/2, conductivity, spoilerlength);
 // spoiler->SetWakePotentials(resWake);


//-----------------------------------------------------
//             Create the WakeProcess
//             and tie it to the wake potential
//-----------------------------------------------------

  SpoilerWakeProcess* collWakeProc = new SpoilerWakeProcess(modes, 1, 100, 3);

  // here we connect the SpoilerWakeProcess  
  // and the TaperedCollimatorPotentials
  // i.e. the SpoilerWakeProcess will only 
  // take care if he finds a wake potential
  // of type TaperedCollimatorPotentials 
  collWake->SetExpectedProcess(collWakeProc);
 
//-----------------------------------------------------
//             Create a wake potential and
//             add it to the cavity element
//-----------------------------------------------------

  TeslaWakePotentials* cavWake = new TeslaWakePotentials;
  cavity->SetWakePotentials(cavWake);

//-----------------------------------------------------
//             Create the WakeProcess and
//             tie it to the wake potential
//-----------------------------------------------------

  WakeFieldProcess* cavWakeProc = new WakeFieldProcess(1);
  cavWakeProc->ApplyImpulseAt(WakeFieldProcess::atExit);
  cavWake->SetExpectedProcess(cavWakeProc);
  
//-----------------------------------------------------
//             Construct Tracker and add
//             the different processes
//-----------------------------------------------------

  ParticleTracker* tracker = new ParticleTracker(model->GetBeamline(),startBunch);

  tracker->AddProcess(collimation);
  tracker->AddProcess(collWakeProc);
  tracker->AddProcess(cavWakeProc);

//-----------------------------------------------------
//             Final bunch parameters
//-----------------------------------------------------
   
  ParticleBunch* finalBunch = tracker->Track(startBunch);

  ofstream final("FinalParameters.dat");
  finalBunch->Output(final);
   
  int i=0;
  double averageyp=0;
  for(ParticleBunch::iterator p=finalBunch->begin();p!=finalBunch->end();p++)
      {averageyp+=(p->yp()-averageyp)/(++i); } 
 
  double xmeanf = finalBunch->GetMoments(0).first;
  double xsigf  = finalBunch->GetMoments(0).second;
  double ymeanf = finalBunch->GetMoments(2).first;
  double ysigf  = finalBunch->GetMoments(2).second;
         
  cout<<endl<<"Final bunch parameters:"<<endl;
  cout<<"Mean x ="<<xmeanf<<' '<<"Sigma x ="<<xsigf<<endl;
  cout<<"Mean y ="<<ymeanf<<' '<<"Sigma y ="<<ysigf<<endl;
  cout<<"yp angle final:"<<finalBunch->GetMoments(3).first<<endl;

  return 0;

}





      

    
  
             
