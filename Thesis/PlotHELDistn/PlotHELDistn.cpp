//Include relevant C++ headers

#include <iostream> // input/output
#include <sstream> // handles string streams
#include <string>
#include <map>
#include <set>
#include <ctime> // used for random seed
#include <sys/stat.h> //to use mkdir

//include relevant MERLIN headers
#include "AcceleratorModel/Apertures/CollimatorAperture.h"

#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchTypes.h"
#include "BeamDynamics/ParticleTracking/HollowELensProcess.h"
#include "BeamDynamics/ParticleTracking/CCFailureProcess.h"
#include "BeamDynamics/ParticleTracking/Integrators/SymplecticIntegrators.h"

#include "Collimators/CollimateProtonProcess.h"
#include "Collimators/ScatteringProcess.h"
#include "Collimators/ScatteringModel.h"
#include "Collimators/CollimatorDatabase.h"
#include "Collimators/MaterialDatabase.h"
#include "Collimators/ApertureConfiguration.h"
#include "Collimators/Dustbin.h"

#include "MADInterface/MADInterface.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"

#include "Random/RandomNG.h"

#include "RingDynamics/Dispersion.h"
#include "RingDynamics/PhaseAdvance.h"
#include "RingDynamics/TransferMatrix.h"
#include "BeamDynamics/ParticleTracking/RingDeltaTProcess.h"

// C++ std namespace, and MERLIN PhysicalUnits namespace

using namespace std;
using namespace PhysicalUnits;

int main(int argc, char* argv[])
{
    int seed = (int)time(NULL);                 // seed for random number generators
    int npart = 1E4;                          	// number of particles to track
    int nturns = 10;
    
    if (argc >=2){	npart = atoi(argv[1]);	}
    if (argc >=3){  seed = atoi(argv[2]);	}
	string bunch_dir, loss_dir, sp_dir, pn_dir, tune_dir, lattice_dir, case_dir, full_case_dir, full_output_dir;
		
    //~ string directory = "/afs/cern.ch/user/h/harafiqu/public/MERLIN";		//lxplus harafiqu
	//~ string directory = "/home/haroon/git/Merlin/HR";						//iiaa1
	string directory = "/home/HR/Downloads/MERLIN";								//M11x	

	string input_dir = "/Thesis/data/PlotHELDistn/";
	
	string output_dir = "/Build/Thesis/outputs/PlotHELDistn/";	

	bool batch = 1;
	if(batch){
		case_dir = "nH/";	
		full_output_dir = (directory+output_dir+case_dir);
		mkdir(full_output_dir.c_str(), S_IRWXU);
	}
	else{		
		string full_output_dir = (directory+output_dir);
		mkdir(full_output_dir.c_str(), S_IRWXU);
	}
	
	bool cut_distn					= 1;
	bool collimation_on 			= 0;
		bool sixtrack 				= 0; 		//use sixtrack like scattering?
		bool dust 					= 0; 		//use dustbin
		if (collimation_on && dust){
			loss_dir = (full_output_dir+"Loss/");
			mkdir(loss_dir.c_str(), S_IRWXU);
		}
		
	bool output_initial_bunch 	= 1;
	bool output_final_bunch 	= 1;
		if (output_initial_bunch || output_final_bunch){
			bunch_dir = (full_output_dir+"Bunch_Distn/");
			mkdir(bunch_dir.c_str(), S_IRWXU);
		}

	bool track 					= 1;	// tracking

	bool hel_on 				= 0; 		//Hollow electron lens process?
	bool LHC_HEL				= 1;
		bool DCon							= 0;
		bool ACon							= 0;
		if(ACon){DCon=0;}
		bool Turnskipon						= 0;
		if(Turnskipon){ACon=0; DCon=0;}
		bool Diffusiveon					= 1;
		if(Diffusiveon){ACon=0; Turnskipon=0; DCon=0;}

	bool every_bunch			= 1;		//output whole bunch every turn in a single file


	/*******************************************************
	** Initialise the random number generator with the seed
	*******************************************************/
    RandomNG::init(seed);
    double beam_energy = 7000.0;

    cout << "npart=" << npart << " nturns=" << nturns << " beam energy = " << beam_energy << endl;

    // Define useful variables
    double beam_charge = 1.1e11;
    double normalized_emittance = 3.5e-6;
    double gamma = beam_energy/PhysicalConstants::ProtonMassMeV/PhysicalUnits::MeV;
	double beta = sqrt(1.0-(1.0/pow(gamma,2)));
	double emittance = normalized_emittance/(gamma*beta);
	
	/********************************
	** ACCELERATORMODEL CONSTRUCTION
	********************************/
	cout << "MADInterface" << endl;
	bool crossing = 0;
	//~ MADInterface* myMADinterface = new MADInterface( directory+input_dir+"twiss_hllhc_b1_thick_HEL.tfs", beam_energy ); 		//HEL
	MADInterface* myMADinterface;
	//~ if(!crossing)
		//~ myMADinterface = new MADInterface( directory+input_dir+"twiss_hllhc_b1_thick_CC.tfs", beam_energy ); 			//noHEL
	//~ else
		myMADinterface = new MADInterface( directory+input_dir+"HLv1.2.0+HEL.tfs", beam_energy ); 	//CrossingOn
	//~ myMADinterface->TreatTypeAsDrift("RFCAVITY");
    //~ myMADinterface->TreatTypeAsDrift("SEXTUPOLE");
    //~ myMADinterface->TreatTypeAsDrift("OCTUPOLE");
	myMADinterface->ConstructApertures(false);
    AcceleratorModel* myAccModel = myMADinterface->ConstructModel();    

	string tcp_element = "TCP.C6L7.B1";    // HORIZONTAL COLLIMATOR (x)
    int tcp_element_number = myAccModel->FindElementLatticePosition(tcp_element.c_str());

    string hel_element = "HEL.B5R4.B1";
	int hel_element_number = myAccModel->FindElementLatticePosition(hel_element.c_str());

	string end_element = "IP1.L1";
	int end_element_number = myAccModel->FindElementLatticePosition(end_element.c_str());
	
	string ip1_element = "IP1";
	int ip1_element_number = myAccModel->FindElementLatticePosition(ip1_element.c_str());

	int start_element_number = tcp_element_number;
		string start_element = tcp_element;
	int plot_element_number = hel_element_number;
		string plot_element = hel_element;
	
    cout << "Found TCP.C6L7 at element number " << tcp_element_number << endl;
    cout << "Plotting distribution at  " << plot_element << " at element number " << plot_element_number << endl;
    cout << "Start element number " << start_element_number << endl;
    cout << "End element number " << end_element_number << endl;

	/*********
	** TWISS
	**********/
	LatticeFunctionTable* myTwiss = new LatticeFunctionTable(myAccModel, beam_energy);
    myTwiss->AddFunction(1,6,3);
    myTwiss->AddFunction(2,6,3);
    myTwiss->AddFunction(3,6,3);
    myTwiss->AddFunction(4,6,3);
    myTwiss->AddFunction(6,6,3);
    myTwiss->AddFunction(0,0,1);
    myTwiss->AddFunction(0,0,2);    
    myTwiss->AddFunction(0,0,3);
    
    double bscale1 = 1e-22;
    while(true){
	cout << "start while(true) to scale bend path length" << endl;
		myTwiss->ScaleBendPathLength(bscale1);
		myTwiss->Calculate();
		if(!std::isnan(myTwiss->Value(1,1,1,0))) {break;}
		bscale1 *= 2;
	}

	/*******************
	** Collimator Setup 
	********************/
	MaterialDatabase* myMaterialDatabase = new MaterialDatabase();
	CollimatorDatabase* collimator_db = new CollimatorDatabase( directory+input_dir+"ColDB_HLv1.2.0.txt", myMaterialDatabase,  true);
   
    collimator_db->MatchBeamEnvelope(true);
    collimator_db->EnableJawAlignmentErrors(false);

    collimator_db->SetJawPositionError(0.0 * nanometer);
    collimator_db->SetJawAngleError(0.0 * microradian);
	//6 sigma HEL Nominal
	//~ collimator_db->SelectImpactFactor("TCP.C6L7.B1", -5.34135E-5);
	//6 sigma HEL HL
	collimator_db->SelectImpactFactor("TCP.C6L7.B1", -5.32764E-5);

	//~ collimator_db->SelectImpactFactor(start_element, 1.0e-6);

    double impact = 1;
    try{ impact = collimator_db->ConfigureCollimators(myAccModel, emittance, emittance, myTwiss);}
    catch(exception& e){ std::cout << "Exception caught: " << e.what() << std::endl; exit(1); }
    if(std::isnan(impact)){ cerr << "Impact is nan" << endl; exit(1); }
    cout << "Impact factor number of sigmas: " << impact << endl;
    delete collimator_db;

	/*************
	** APERTURES
	**************/
	ApertureConfiguration* myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_HLv1.2.0.tfs",1);      
    
    myApertureConfiguration->ConfigureElementApertures(myAccModel);
    delete myApertureConfiguration;

	//CHECK FOR COLLIMATOR APERTURES	
	vector<Collimator*> TCP;
	int siz = myAccModel->ExtractTypedElements(TCP, start_element);
	
	cout << "\n\t Found " << TCP.size() << " Collimators when extracting" << endl;

	Aperture *ap = (TCP[0])->GetAperture();
	if(!ap)
		{cout << "Could not get tcp ap" << endl;	abort();}
	else{cout << "TCP aperture type = " << ap->GetApertureType() << endl;}

	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
	if(!CollimatorJaw)
		{cout << "Could not cast" << endl;	abort();}

	/********
	** BEAM 
	*********/
	Dispersion* myDispersion = new Dispersion(myAccModel, beam_energy);
    myDispersion->FindDispersion(start_element_number);

    BeamData mybeam;

    mybeam.charge = beam_charge/npart;
    mybeam.p0 = beam_energy;
    mybeam.beta_x = myTwiss->Value(1,1,1,start_element_number)*meter;
    mybeam.beta_y = myTwiss->Value(3,3,2,start_element_number)*meter;
    mybeam.alpha_x = -myTwiss->Value(1,2,1,start_element_number);
    mybeam.alpha_y = -myTwiss->Value(3,4,2,start_element_number);

    // Dispersion
    mybeam.Dx=myDispersion->Dx;
    mybeam.Dy=myDispersion->Dy;
    mybeam.Dxp=myDispersion->Dxp;
    mybeam.Dyp=myDispersion->Dyp;

    //~ mybeam.emit_x = impact * impact * emittance * meter;
    //~ impact =1;
    //~ mybeam.emit_y = impact * impact * emittance * meter;

    //~ mybeam.emit_x = sqrt( (em_t*myTwiss->Value(1,1,1,start_element_number)*meter)/(beta*gamma) );
    //~ mybeam.emit_y = sqrt( (em_t*myTwiss->Value(3,3,2,start_element_number)*meter)/(beta*gamma) );
    mybeam.emit_y = emittance * meter;
    mybeam.emit_x = emittance * meter;
    
    //LHC beam
    mybeam.sig_z = 7.55E-2;
    mybeam.sig_dp = 1.13E-4;

    //Beam centroid
    mybeam.x0=myTwiss->Value(1,0,0,start_element_number);
    mybeam.xp0=myTwiss->Value(2,0,0,start_element_number);
    mybeam.y0=myTwiss->Value(3,0,0,start_element_number);
    mybeam.yp0=myTwiss->Value(4,0,0,start_element_number);
    mybeam.ct0=myTwiss->Value(5,0,0,start_element_number);

    // X-Y coupling
    mybeam.c_xy=0.0;
    mybeam.c_xyp=0.0;
    mybeam.c_xpy=0.0;
    mybeam.c_xpyp=0.0;

    delete myDispersion;

	/********
	** BUNCH 
	*********/
	ProtonBunch* myBunch;
    int node_particles = npart;
    ParticleBunchConstructor* myBunchCtor;
	
    // horizontalHaloDistribution1 is a halo in xx' plane, zero in yy'
    // horizontalHaloDistribution2 is a halo in xx' plane, gaussian in yy'
	//~ myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, CCDistn);
	//~ myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, CCDistn2);
	//~ myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, RFDistn);
	myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, horizontalHaloDistribution1);

	if(cut_distn){ 
			double h_offset = myTwiss->Value(1,0,0,start_element_number);
			double JawPosition = (CollimatorJaw->GetFullWidth() / 2.0);
			cout << "h_offset: " << h_offset << endl;	
			cout << "Jaw position: " << JawPosition << endl;

			HorizontalHaloParticleBunchFilter* hFilter = new HorizontalHaloParticleBunchFilter();
			double tcpsig = 267.5E-6;
			hFilter->SetHorizontalLimit(4*tcpsig);
			hFilter->SetHorizontalOrbit(h_offset);
			
			myBunchCtor->SetFilter(hFilter); 
	}

    myBunch = myBunchCtor->ConstructParticleBunch<ProtonBunch>();
    delete myBunchCtor;

    myBunch->SetMacroParticleCharge(mybeam.charge);
		
		/***********************
		** Output Initial Bunch
		***********************/   
		ostringstream bunch_output_file;
		bunch_output_file << (bunch_dir + "initial_bunch.txt");

		ofstream* bunch_output = new ofstream(bunch_output_file.str().c_str());
		myBunch->Output(*bunch_output);
		delete bunch_output;

	/******************
	** ParticleTracker
	******************/
	ParticleTracker* myPreParticleTracker;
	ParticleTracker* myPrePreParticleTracker;
	ParticleTracker* myParticleTracker;

	AcceleratorModel::Beamline preprebeamline = myAccModel->GetBeamline(tcp_element_number, end_element_number);
	AcceleratorModel::Beamline prebeamline = myAccModel->GetBeamline(ip1_element_number, hel_element_number-1);
	AcceleratorModel::RingIterator beamline = myAccModel->GetRing(hel_element_number);
    
	myPrePreParticleTracker = new ParticleTracker(preprebeamline, myBunch);
	myPreParticleTracker = new ParticleTracker(prebeamline, myBunch);
	myParticleTracker = new ParticleTracker(beamline, myBunch);
	
	myPrePreParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());		
	myPreParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());		
	myParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());	
	//~ myPrePreParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());	
	//~ myPreParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());		
	//~ myParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());

	/***************
	** COLLIMATION
	****************/
	LossMapDustbin* myDustbin = new LossMapDustbin;
	if(collimation_on){
		CollimateProtonProcess* myCollimateProcess =new CollimateProtonProcess(2, 4, NULL);
	
		myCollimateProcess->SetDustbin(myDustbin);       
			//~ myCollimateProcess->ScatterAtCollimator(false);
			myCollimateProcess->ScatterAtCollimator(true);
		ScatteringModel* myScatter = new ScatteringModel;
	  
		if(sixtrack){
			myScatter->SetScatterType(0);
		}
		else{
			myScatter->SetScatterType(4);
		}

		myCollimateProcess->SetScatteringModel(myScatter);
		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);

		myPrePreParticleTracker->AddProcess(myCollimateProcess);
		myPreParticleTracker->AddProcess(myCollimateProcess);
		myParticleTracker->AddProcess(myCollimateProcess);
	}			
	
	/***************
	** HEL PROCESS
	****************/
	if(hel_on){
		
		// HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double lengh_e);
		HollowELensProcess* myHELProcess;
			
		if(LHC_HEL){	// LHC: 3m, 10KeV, 5A
			myHELProcess = new HollowELensProcess(3, 1, 5, 0.195, 2.334948339E4, 3.0);
		}
		else{			//Tevatron: 2m, 5KeV, 1.2A
			myHELProcess = new HollowELensProcess(2, 1, 1.2, 0.138874007, 2.334948339E4, 2.0);
		}
				
		myHELProcess->SetRadialProfile();
		//~ myHELProcess->SetPerfectProfile();
		
		if(LHC_HEL){
			myHELProcess->SetRadiiSigma(4, 8, myAccModel, emittance, emittance, myTwiss);
		}
		else{//Tevatron
			myHELProcess->SetRadiiSigma(4, 6.8, myAccModel, emittance, emittance, myTwiss);
		}		
		
		if(ACon){
			//Set AC variables
			//HollowELensProcess::SetAC (double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi) 
			//H20
			myHELProcess->SetAC(0.31, .002, 5E-5, 1E3, 2.);
			myHELProcess->SetOpMode(AC);
		}
		else if(DCon){
			myHELProcess->SetOpMode(DC);
		}
		else if(Diffusiveon){	
			myHELProcess->SetOpMode(Diffusive);
		}
		else if(Turnskipon){	
			myHELProcess->SetTurnskip(2);
			myHELProcess->SetOpMode(Turnskip);
		}
		//~ myPreParticleTracker->AddProcess(myHELProcess);	
		myParticleTracker->AddProcess(myHELProcess);	
	}
	
	// Output bunch every turn in one file
	ostringstream bo_file;
	bo_file << bunch_dir << "total_bunch.txt";
	
	//truncate (clear) the file first to prevent appending to last run
	ofstream* boclean = new ofstream(bo_file.str().c_str(), ios::trunc);
	ofstream* bo = new ofstream(bo_file.str().c_str(), ios::app);
	
	/***************
	** TRACKING RUN
	****************/
	
	
	myPrePreParticleTracker->Track(myBunch);
	myPreParticleTracker->Track(myBunch);

	for (int turn=1; turn<=nturns; turn++)
	{
		cout << "Start tracking loop " << turn << endl;
		
		myParticleTracker->Track(myBunch);
		
		cout << "PreTracked " << turn << " turns" << endl;
		if(every_bunch){myBunch->Output(*bo);}
		

		if( myBunch->size() <= 1 ) break;
	}
	/*****************
	** DUSTBIN LOSSES
	******************/
	if(dust){
		ostringstream dustbin_file;
		dustbin_file << (loss_dir+"HELdustbin_losses.txt");	
		ofstream* dustbin_output = new ofstream(dustbin_file.str().c_str());
		if(!dustbin_output->good())    {
			std::cerr << "Could not open dustbin loss file" << std::endl;
			exit(EXIT_FAILURE);
		} 
		
		myDustbin->Finalise(); 
		myDustbin->Output(dustbin_output); 
	}			
	/***************
	** FINAL DISTN
	****************/
	if(output_final_bunch){
		ostringstream bunch_output_file2;
		bunch_output_file2 << bunch_dir << "final_bunch.txt";

		ofstream* bunch_output2 = new ofstream(bunch_output_file2.str().c_str());
		myBunch->Output(*bunch_output2);
		delete bunch_output2;
	 }
   
    cout << "npart: " << npart << endl;
    cout << "left: " << myBunch->size() << endl;
    cout << "absorbed: " << npart - myBunch->size() << endl;

    delete myMaterialDatabase;
    delete myBunch;
    delete myTwiss;
    delete myAccModel;
    delete myMADinterface;	

    return 0;
}
