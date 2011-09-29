//------------------------------------------------------------
//
// an example how to use ROOT with Merlin
//
// We construct an ILC like linac and
// use the TrackingOutputROOT class 
// to write all kinds of variables and parameters
//
// D.K. 14.9.2009
//

#include "merlin_config.h"
#include "Random/RandomNG.h"
#include "MADInterface/XTFFInterface.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/BeamData.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"

#include "NumericalUtils/NumericalConstants.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/PhysicalUnits.h"

#include <iostream> 

#include "TrackingOutputROOT.h"
#include "TeslaWakePotential.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname);

int main() {

//-------------------------- an example MERLIN simulation -----------------------------------

	// Initialise the random number generator
	RandomNG::init(1);

	// Construct model
	pair<AcceleratorModel*,BeamData*> mb = ConstructModel("ilc_linac_15_250.xtff");  

	AcceleratorModel*          model = mb.first;
	BeamData*                  beam  = mb.second;
	AcceleratorModel::Beamline bline = model->GetBeamline();


//-------------------------- here we demonstrate the ROOT interface -------------------------

	// output some bunch/lattice parameters along the track
	//
	TrackingOutputROOT* trackout = new TrackingOutputROOT();
	//
	// where we want to get the output
	//
	trackout->output_all     = false;
	trackout->output_final   = true; 
	trackout->output_initial = true;
	trackout->AddIdentifier("Quadrupole.*");
//	trackout->AddIdentifier("BPM.*");
//	trackout->AddIdentifier("YCor.*");
//	trackout->AddIdentifier("TWRFStructure.*");
//	trackout->AddIdentifierAllMags();
	//
	// what kind of data
	trackout->Output.base=true;         // component z, reference momentum p0 (bunch E=p0(1+dp))
	trackout->Output.names=true;        // type, name
	trackout->Output.bunchFirst=true;   // means      - x,xp,y,yp,ct,dp   : <x>,<x'> etc.
	trackout->Output.bunchSecond=true;  // sigmas     - s_x,s_xp,s_y,s_yp,s_ct,s_dp   : sigma(x) etc.
	trackout->Output.emittance=true;    // emittances - gex,gey,gexc,geyc : c energy corr.
	trackout->Output.Twiss=true;        // twiss parameters - ax,bx,dx,dxp,ay,by,dy,dyp
	trackout->Output.BPMCor=true;       // BPM readings and Corrector settings - BPM_x BPM_y, XCor, YCor
	trackout->Output.Cav=true;          // Cavity voltage and phase
//	trackout->Output.AllFalse();      // Switch off all
	
	// dump the bunch at id
	//
	// bunch tree name = track tree name + name part of id 
	//                                 ( + counter if more than one)
	trackout->DumpBunchAt("INITIAL");
	trackout->DumpBunchAt("FINAL");
//	trackout->DumpBunchAt("BPM.*");  // at each! BPM, large output

	// A tracker
	ParticleTracker tracker(bline);
	//
	// with wakefield
 	WakeFieldProcess* wakep = new WakeFieldProcess(1);
 	wakep->ApplyImpulseAt(WakeFieldProcess::atExit);
 	tracker.AddProcess(wakep);

	// here we plug in the ROOT output !
	tracker.SetOutput(trackout);

	// We run 3 different bunches through the tracker
	// Each bunch gets its own TrackingTree and BunchTree with specific names 
	for(int i=0;i<3;i++){ 	

		cout<<endl<<"Bunch "<<i<<endl;
		
		// we restrict the ROOT output for the last bunch
		if(i==2){
			trackout->Output.AllFalse();
			trackout->Output.names=true;
			trackout->Output.base=true;
			trackout->Output.emittance=true;
			// no bunch output
			trackout->DoBunch(false);
		}
		
		// define a unique name
		stringstream treename;
		treename<<"mytree_"<<i;
		trackout->NewTree(treename.str());
				
		// Construct bunch according to beam parameter
		ParticleBunch* bunch = ParticleBunchConstructor(*beam,1000).ConstructParticleBunch();
		
		cout<<"Initial energy "<<bunch->GetReferenceMomentum()<<" GeV"<<endl;
		tracker.Track(bunch);
		cout<<"Final energy "<<bunch->GetReferenceMomentum()<<" GeV"<<endl;

		delete bunch;
	
	}	
	delete model;
	delete beam;

	delete trackout; // necessary for output

	return 0;
}

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname) {

	double qt = 2.0e+10; // single-bunch charge

	ofstream logs("construction.log");
	XTFFInterface mc(fname,qt,&logs);
	mc.ConstructGirders(true);

	pair<AcceleratorModel*,BeamData*> mb = mc.Parse();
	BeamData* beam0 = mb.second;
	AcceleratorModel* model = mb.first;
	
	// The following quantities are not
	// given by the XTFF file
	
	double gamma = beam0->p0/MeV/ElectronMassMeV;
	beam0->emit_x = 8.0e-06/gamma;
	beam0->emit_y = 0.02e-06/gamma;
	beam0->charge = qt==0 ? 1.0 : qt;
	beam0->sig_dp = 0.0107;
	beam0->sig_z  = 300.0e-06;

	// include wakefields
	vector<TWRFStructure*> cavities;
	int nc = model->ExtractTypedElements(cavities,"CAV");
	TeslaWakePotentials* wake = new TeslaWakePotentials;
        for(vector<TWRFStructure*>::iterator c = cavities.begin(); c!=cavities.end(); c++) {
                (*c)->SetWakePotentials(wake);
        }
        cout<<nc<<" cavities with Wakefields."<<endl; 

	return mb;

}


