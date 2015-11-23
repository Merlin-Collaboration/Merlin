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

// C++ std namespace, and MERLIN PhysicalUnits namespace

using namespace std;
using namespace PhysicalUnits;

// Main function, this executable can be run with the arguments number_of_particles seed

//e.g. for 1000 particles and a seed of 356: ./test 1000 356

int main(int argc, char* argv[])
{
    int seed = (int)time(NULL);                 // seed for random number generators
    int npart = 51;                            // number of particles to track
    int nturns = 1E4;                           // number of turns to track
 
    if (argc >=2){
        npart = atoi(argv[1]);
    }

    if (argc >=3){
        seed = atoi(argv[2]);
    }

    RandomNG::init(seed);
    double beam_energy = 7000.0;
    cout << "npart=" << npart << " nturns=" << nturns << " beam energy = " << beam_energy << endl;

    // Define useful variables
    double beam_charge = 1.1e11;
    double normalized_emittance = 3.5e-6;
    double gamma = beam_energy/PhysicalConstants::ProtonMassMeV/PhysicalUnits::MeV;
	double beta = sqrt(1.0-(1.0/pow(gamma,2)));
	double emittance = normalized_emittance/(gamma*beta);
	
	//~ string directory = "/afs/cern.ch/user/h/harafiqu/public/MERLIN";	//lxplus harafiqu
	//~ string directory = "/home/haroon/git/Merlin";				//iiaa1
	string directory = "/home/HR/Downloads/MERLIN";					//M11x	
	//~ string directory = "/afs/cern.ch/user/a/avalloni/private/Merlin_all";	//lxplus avalloni
	
	string pn_dir, case_dir, bunch_dir;	
	string input_dir = "/Thesis/data/HELIntegration/";	
	string output_dir = "/Build/Thesis/outputs/HELIntegration/";
	
	bool batch = 1;
	if(batch){
		//~ case_dir = "DC_PhaseEllipse_-90m/";
		case_dir = "Poincare_-30m/";
	}
	
	string full_output_dir = (directory+output_dir);
	mkdir(full_output_dir.c_str(), S_IRWXU);
	
	if(batch){
		full_output_dir = (directory+output_dir+case_dir);
		mkdir(full_output_dir.c_str(), S_IRWXU);
	}
	bool output_turn_bunch		= 1;
		if(output_turn_bunch){
			pn_dir = (full_output_dir+"ParticleNo/");
			mkdir(pn_dir.c_str(), S_IRWXU);
	}		
	bool every_bunch			= 1;		//output whole bunch every turn in a single file
	bool output_initial_bunch 	= 1;
	bool output_final_bunch 	= 1;
		if (output_initial_bunch || output_final_bunch){
			bunch_dir = (full_output_dir+"Bunch_Distn/");
			mkdir(bunch_dir.c_str(), S_IRWXU);
		}		
	bool output_fluka_database 	= 1;
	
	bool hel_on 				= 1; 		//Hollow electron lens process?
		bool DCon							= 1;
		bool ACon							= 0;
		if(ACon){DCon=0;}
		bool Turnskipon						= 0;
		if(Turnskipon){ACon=0; DCon=0;}
		bool Diffusiveon					= 0;
		if(Diffusiveon){ACon=0; Turnskipon=0; DCon=0;}
		
	bool collimation_on 		= 0;
	bool cut_distn				= 0;
	bool round_beams			= 1;
	
	//REMEMBER TO CHANGE DISTRIBUTION SIGMA
	// note that this gives the correct phase advance as we don't use m.apply()
	bool start_at_ip1			= 1;	// True: 3 trackers: IP1->HEL, HEL->TCP, TCP->IP1 
										// False: 3 trackers: TCP->IP1, IP1->HEL, HEL->TCP
										
	bool cleaning				= 0;
		if(cleaning){
			collimation_on		= 1;
			every_bunch			= 0;
			output_turn_bunch	= 1;
			start_at_ip1		= 0;	
			cut_distn			= 1;	
			output_initial_bunch= 1;
			output_final_bunch	= 1;
		}
		
///////////////////////////////////
// ACCELERATORMODEL CONSTRUCTION //
///////////////////////////////////
	cout << "MADInterface" << endl;
	MADInterface* myMADinterface;
	
    //~ MADInterface* myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_RF.tfs", beam_energy );
    //~ MADInterface* myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_original.tfs", beam_energy );
    if(round_beams)
		myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_RF_-30mHEL.tfs", beam_energy );	//new HL v1.2
		//~ myMADinterface = new MADInterface( directory+input_dir+"HLv1.2.0_C+S_RF_-30mHEL.tfs", beam_energy );		//old HL v1.2
    else
		myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_RF_-90mHEL.tfs", beam_energy );

    // As we are only tracking for 200 turns we can choose to ignore the accelerating cavities
    // To do this we use the TreatTypeAsDrift() function, which takes an element type string as an argument, this can be done for any element
    //~ myMADinterface->TreatTypeAsDrift("RFCAVITY");
    //~ myMADinterface->TreatTypeAsDrift("SEXTUPOLE");
    //~ myMADinterface->TreatTypeAsDrift("OCTUPOLE");

    myMADinterface->ConstructApertures(false);

    AcceleratorModel* myAccModel = myMADinterface->ConstructModel();   

///////////
// TWISS //
///////////

	int hel_element_number = 0;
	string hel_element;
	// Primary collimator
	string tcp_element = "TCP.C6L7.B1";    // HORIZONTAL COLLIMATOR (x)
    int tcp_element_number = myAccModel->FindElementLatticePosition(tcp_element.c_str()); 
    
    // Hollow electron lens
    if(round_beams){
		hel_element = "HEL_-30.B5L4.B1";
		hel_element_number = myAccModel->FindElementLatticePosition(hel_element.c_str());
	}
    else{
		hel_element = "HEL_-90.A5L4.B1";
		hel_element_number = myAccModel->FindElementLatticePosition(hel_element.c_str());
	}
    // END of Lattice
    string end_element = "IP1.L1";
	int end_element_number = myAccModel->FindElementLatticePosition(end_element.c_str());
	// START of Lattice
    string ip1_element = "IP1";
	int ip1_element_number = myAccModel->FindElementLatticePosition(ip1_element.c_str());	
    
    int start_element_number;
    if(start_at_ip1)
		start_element_number = ip1_element_number;
	else
		start_element_number = tcp_element_number;

    cout << "Found start element IP1 at element number " << start_element_number << endl;
    cout << "Found start element HEL at element number " << hel_element_number << endl;
    cout << "Found start element TCP.C6L7 at element number " << tcp_element_number << endl;
    cout << "Found start element END at element number " << end_element_number << endl;
    cout << "Found start element START at element number " << start_element_number << endl;

  LatticeFunctionTable* myTwiss = new LatticeFunctionTable(myAccModel, beam_energy);
    myTwiss->AddFunction(1,6,3);
    myTwiss->AddFunction(2,6,3);
    myTwiss->AddFunction(3,6,3);
    myTwiss->AddFunction(4,6,3);
    myTwiss->AddFunction(6,6,3);

    double bscale1 = 1e-22;
    
  
	while(true)
	{
	cout << "start while(true) to scale bend path length" << endl;
		myTwiss->ScaleBendPathLength(bscale1);
		myTwiss->Calculate();

		if(!std::isnan(myTwiss->Value(1,1,1,0))) {break;}
		bscale1 *= 2;
	}
	
	Dispersion* myDispersion = new Dispersion(myAccModel, beam_energy);
    myDispersion->FindDispersion(start_element_number);

	string lattice_dir = (full_output_dir+"LatticeFunctions/");
	mkdir(lattice_dir.c_str(), S_IRWXU);
	
	ostringstream twiss_output_file; 
    twiss_output_file << (lattice_dir+"LatticeFunctions.dat");
    ofstream twiss_output(twiss_output_file.str().c_str());
	if(!twiss_output.good())
	{
		std::cerr << "Could not open twiss output file" << std::endl;
		exit(EXIT_FAILURE);
	} 
	myTwiss->PrintTable(twiss_output);
	
	// Doesn't work (yet)
	ostringstream disp_output_file; 
    disp_output_file << (lattice_dir+"Dispersion.dat");
    ofstream* disp_output = new ofstream(disp_output_file.str().c_str());
	if(!disp_output->good())
	{
		std::cerr << "Could not open dispersion output file" << std::endl;
		exit(EXIT_FAILURE);
	} 
	//~ myDispersion->OutputDispersion(disp_output);
	myDispersion->FindRMSDispersion(disp_output);
	
///////////////////////
// Collimator set up //
///////////////////////
	cout << "Collimator Setup" << endl;   
   
    MaterialDatabase* myMaterialDatabase = new MaterialDatabase();
    CollimatorDatabase* collimator_db = new CollimatorDatabase( directory+input_dir+"HL_v1.2.1_collDB.txt", myMaterialDatabase,  true);
   
    collimator_db->MatchBeamEnvelope(true);
    collimator_db->EnableJawAlignmentErrors(false);
    collimator_db->SetJawPositionError(0.0 * nanometer);
    collimator_db->SetJawAngleError(0.0 * microradian);
    
    // HLv1.2  -0.1 sigma = -2.73539E-5
    // HLv1.2  -0.2 sigma = -5.328E-5
    // HLv1.2  -0.3 sigma = -7.992E-5
    // HLv1.2  -1.2 sigma = -3.196E-4
    // HLv1.2  -0.7 sigma = -1.8648E-4
    
    //~ collimator_db->SelectImpactFactor(tcp_element, -2.66313E-5);    
    collimator_db->SelectImpactFactor(tcp_element, -(10*2.66313E-5));    
	//~ collimator_db->SelectImpactFactor(tcp_element, 1.0e-6);

    double impact = 6;
    
    try{
        impact = collimator_db->ConfigureCollimators(myAccModel, emittance, emittance, myTwiss);
    }
    catch(exception& e){
        std::cout << "Exception caught: " << e.what() << std::endl;
        exit(1);
    }
    if(std::isnan(impact)){
        cerr << "Impact is nan" << endl;
        exit(1);
    }
    cout << "Impact factor number of sigmas: " << impact << endl;
    
    if(output_fluka_database){
		ostringstream fd_output_file;
		fd_output_file << (full_output_dir+"fluka_database.txt");

		ofstream* fd_output = new ofstream(fd_output_file.str().c_str());
		collimator_db->OutputFlukaDatabase(fd_output);
		delete fd_output;
	}    
    
    delete collimator_db;
    
//CHECK FOR COLLIMATOR APERTURES	
	vector<Collimator*> TCP;
	int siz = myAccModel->ExtractTypedElements(TCP, tcp_element);

	cout << "\n\t Found " << TCP.size() << " Collimators when extracting" << endl;

	Aperture *ap = (TCP[0])->GetAperture();
	if(!ap){cout << "Could not get tcp ap" << endl;	abort();}
	else{cout << "TCP aperture type = " << ap->GetApertureType() << endl;}

	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
	if(!CollimatorJaw){cout << "Could not cast" << endl;	abort();}
	
////////////////////////////
// Aperture Configuration //
////////////////////////////

	ApertureConfiguration* myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"HL_v1.2.1_Aperture.tfs",1);      
    
    myApertureConfiguration->ConfigureElementApertures(myAccModel);
    delete myApertureConfiguration;

///////////////////
// BEAM SETTINGS //
///////////////////

    BeamData mybeam;

    // Default values for all members of BeamData are 0.0
    // Particles are treated as macro particles - this has no bearing on collimation
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

    // We set the beam emittance such that the bunch (created from this BeamData object later) will impact upon the primary collimator
    //~ mybeam.emit_x = impact * impact * emittance * meter;
    if(start_at_ip1){
		mybeam.emit_x = emittance * meter;
		mybeam.emit_y = emittance * meter;
	}
    else{
		//~ mybeam.emit_x = emittance * meter;
		//~ mybeam.emit_y = emittance * meter;
		mybeam.emit_x = impact * impact * emittance * meter;
		mybeam.emit_y = impact * impact * emittance * meter;
	}
    impact =1;
    //~ mybeam.emit_y = impact * impact * emittance * meter;
    
    mybeam.sig_z = 0.0;

    //Beam centroid
    mybeam.x0=myTwiss->Value(1,0,0,start_element_number);
    mybeam.xp0=myTwiss->Value(2,0,0,start_element_number);
    mybeam.y0=myTwiss->Value(3,0,0,start_element_number);
    mybeam.yp0=myTwiss->Value(4,0,0,start_element_number);
    mybeam.ct0=myTwiss->Value(5,0,0,start_element_number);

    mybeam.sig_dp = 0.0;

    // X-Y coupling
    mybeam.c_xy=0.0;
    mybeam.c_xyp=0.0;
    mybeam.c_xpy=0.0;
    mybeam.c_xpyp=0.0;

    delete myDispersion;

///////////
// BUNCH //
///////////

    ProtonBunch* myBunch;
    int node_particles = npart;

    // horizontalHaloDistribution1 is a halo in xx' plane, zero in yy'
    // horizontalHaloDistribution2 is a halo in xx' plane, gaussian in yy'
    ParticleBunchConstructor* myBunchCtor;
    if(cleaning)
		myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, horizontalHaloDistribution2);
    else
    	myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, tuneTestDistribution);
    
	if(collimation_on && cut_distn){ 
		double h_offset = myTwiss->Value(1,0,0,start_element_number);
		double JawPosition = (CollimatorJaw->GetFullWidth() / 2.0);
		cout << "h_offset: " << h_offset << endl;	
		cout << "Jaw position: " << JawPosition << endl;

		HorizontalHaloParticleBunchFilter* hFilter = new HorizontalHaloParticleBunchFilter();
		//~ double tcpsig = 0.000273539;
		double tcpsig = 0.000266313;
		hFilter->SetHorizontalLimit(4*tcpsig);
		hFilter->SetHorizontalOrbit(h_offset);

		myBunchCtor->SetFilter(hFilter); 
	}
	
	myBunch = myBunchCtor->ConstructParticleBunch<ProtonBunch>();
    delete myBunchCtor;

    myBunch->SetMacroParticleCharge(mybeam.charge);
    
   if(output_initial_bunch){
		ostringstream bunch_output_file;
		bunch_output_file << (bunch_dir + "initial_bunch.txt");
		ofstream* bunch_output = new ofstream(bunch_output_file.str().c_str());
		if(!bunch_output->good())
		{
			std::cerr << "Could not open dustbin loss file" << std::endl;
			exit(EXIT_FAILURE);
		}   
		myBunch->Output(*bunch_output);			
		delete bunch_output;
	}

/////////////////////
// ParticleTracker //
/////////////////////

	//~ AcceleratorModel::RingIterator beamline = myAccModel->GetRing(start_element_number);
    //~ ParticleTracker* myParticleTracker = new ParticleTracker(beamline, myBunch);
    
    ParticleTracker* myParticleTracker1;
    ParticleTracker* myParticleTracker2;
    ParticleTracker* myParticleTracker3;
 /*    
    if(start_at_ip1){
		//~ AcceleratorModel::Beamline beamline1 = myAccModel->GetBeamline(start_element_number, hel_element_number-1);
		AcceleratorModel::Beamline beamline1 = myAccModel->GetBeamline(start_element_number, hel_element_number-1);
		AcceleratorModel::Beamline beamline2 = myAccModel->GetBeamline(hel_element_number, end_element_number);
		
		myParticleTracker1 = new ParticleTracker(beamline1, myBunch);
		myParticleTracker2 = new ParticleTracker(beamline2, myBunch);	
		
		myParticleTracker1->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());	
		myParticleTracker2->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
		//~ myParticleTracker1->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());	
		//~ myParticleTracker2->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());		
	}	
*/
	if(start_at_ip1){
		AcceleratorModel::Beamline beamline1 = myAccModel->GetBeamline(start_element_number, hel_element_number-1);
		AcceleratorModel::Beamline beamline2 = myAccModel->GetBeamline(hel_element_number, tcp_element_number-1);
		AcceleratorModel::Beamline beamline3 = myAccModel->GetBeamline(tcp_element_number, end_element_number);
		
		myParticleTracker1 = new ParticleTracker(beamline1, myBunch);
		myParticleTracker2 = new ParticleTracker(beamline2, myBunch);
		myParticleTracker3 = new ParticleTracker(beamline3, myBunch);
				
		myParticleTracker1->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());	
		myParticleTracker2->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
		myParticleTracker3->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
		//~ myParticleTracker1->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());	
		//~ myParticleTracker2->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
		//~ myParticleTracker3->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());		
	}
	else{
		AcceleratorModel::Beamline beamline1 = myAccModel->GetBeamline(tcp_element_number, end_element_number);
		AcceleratorModel::Beamline beamline2 = myAccModel->GetBeamline(ip1_element_number, hel_element_number-1);
		AcceleratorModel::Beamline beamline3 = myAccModel->GetBeamline(hel_element_number, tcp_element_number-1);
		
		myParticleTracker1 = new ParticleTracker(beamline1, myBunch);
		myParticleTracker2 = new ParticleTracker(beamline2, myBunch);
		myParticleTracker3 = new ParticleTracker(beamline3, myBunch);
				
		myParticleTracker1->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());	
		myParticleTracker2->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
		myParticleTracker3->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
		//~ myParticleTracker1->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());	
		//~ myParticleTracker2->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
		//~ myParticleTracker3->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
	}

/////////////////////////
// Collimation Process //
/////////////////////////
	
	LossMapDustbin* myDustbin = new LossMapDustbin;

	if(collimation_on){
		cout << "Collimation on" << endl;
		CollimateProtonProcess* myCollimateProcess = new CollimateProtonProcess(2, 4);
			
		myCollimateProcess->SetDustbin(myDustbin);       

		myCollimateProcess->ScatterAtCollimator(true);
	   
		ScatteringModel* myScatter = new ScatteringModel;

		// 0: ST,    1: ST + Adv. Ionisation,    2: ST + Adv. Elastic,    3: ST + Adv. SD,     4: MERLIN
		bool use_sixtrack_like_scattering = 0;
		if(use_sixtrack_like_scattering){
			myScatter->SetScatterType(0);
		}
		else{
			myScatter->SetScatterType(4);
		}

		myCollimateProcess->SetScatteringModel(myScatter);

		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);
		
		if(start_at_ip1){
			myParticleTracker1->AddProcess(myCollimateProcess);
			myParticleTracker2->AddProcess(myCollimateProcess);
			myParticleTracker3->AddProcess(myCollimateProcess);
			cout << "Collimation added to tracker" << endl;
		}
		else{
			myParticleTracker1->AddProcess(myCollimateProcess);
			myParticleTracker2->AddProcess(myCollimateProcess);
			myParticleTracker3->AddProcess(myCollimateProcess);
			cout << "Collimation added to tracker" << endl;
		}
	}
    
/////////////////
// HEL Process //
/////////////////	
	if(hel_on){
		cout << "HEL on" << endl;
				
		// HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e);
		HollowELensProcess* myHELProcess;
			
		// LHC: 3m, 10KeV, 5A
		myHELProcess = new HollowELensProcess(3, 1, 5, 0.195, 2.334948339E4, 3.0);		
				
		myHELProcess->SetRadialProfile();
		//~ myHELProcess->SetPerfectProfile();
		
		myHELProcess->SetRadiiSigma(4, 8, myAccModel, emittance, emittance, myTwiss);
			
		
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
		
		if(start_at_ip1){
			myParticleTracker2->AddProcess(myHELProcess);	
			cout << "HEL set" << endl;		
		}
		else{			
			myParticleTracker3->AddProcess(myHELProcess);	
			cout << "HEL set" << endl;		
		}
	}

	/*****************************
	 *  Other Output Files
	 ****************************/
 
	ostringstream particle_no_file;
	particle_no_file << pn_dir<< "No.txt";
	ofstream* particle_no_output = new ofstream(particle_no_file.str().c_str());	
	if(!particle_no_output->good())
	{
		std::cerr << "Could not open particle_no_output file" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// Output bunch every turn @HEL in one file
	ostringstream bo_file;
	bo_file << bunch_dir << "HEL_bunch.txt";
	
	//truncate (clear) the file first to prevent appending to last run
	ofstream* boclean = new ofstream(bo_file.str().c_str(), ios::trunc);
	ofstream* bo = new ofstream(bo_file.str().c_str(), ios::app);	
	if(!bo->good())
	{
		std::cerr << "Could not open every bunch HEL output file" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// Output bunch every turn @TCP in one file
	ostringstream bot_file;
	bot_file << bunch_dir << "TCP_bunch.txt";
	
	//truncate (clear) the file first to prevent appending to last run
	ofstream* botclean = new ofstream(bot_file.str().c_str(), ios::trunc);
	ofstream* bot = new ofstream(bot_file.str().c_str(), ios::app);	
	if(!bot->good())
	{
		std::cerr << "Could not open every bunch TCP output file" << std::endl;
		exit(EXIT_FAILURE);
	}
		
//////////////////
// TRACKING RUN //
//////////////////
	
	if(output_turn_bunch)
		(*particle_no_output) << "0\t" << myBunch->size() << endl;
    
    for (int turn=1; turn<=nturns; turn++)
    {
        cout << "Turn " << turn <<"\tParticle number: " << myBunch->size() << endl;

		if(start_at_ip1){	myParticleTracker1->Track(myBunch);}
		else{				myParticleTracker1->Track(myBunch);
							myParticleTracker2->Track(myBunch);}
							
		if(every_bunch){myBunch->Output(*bo);} //Split the tracker to output at HEL
        
        if(start_at_ip1){	myParticleTracker2->Track(myBunch);}
		else{				myParticleTracker3->Track(myBunch);}	
		
		if(every_bunch){myBunch->Output(*bot);} //Split the tracker to output at TCP
		
		if(start_at_ip1){	myParticleTracker3->Track(myBunch);}	
		
		if(output_turn_bunch){(*particle_no_output) << turn <<"\t" << myBunch->size() << endl;}
			
        if( myBunch->size() <= 1 ) break;
    }
    
////////////////////////////
// OUTPUT DUSTBIN LOSSES ///
////////////////////////////
	if(collimation_on){
		ostringstream dustbin_file;
		dustbin_file << (full_output_dir+"dustbin_losses.txt");	
		ofstream* dustbin_output = new ofstream(dustbin_file.str().c_str());	
		if(!dustbin_output->good())    {
			std::cerr << "Could not open dustbin loss file" << std::endl;
			exit(EXIT_FAILURE);
		}   
		myDustbin->Finalise(); 
		myDustbin->Output(dustbin_output); 
	}
	
//////////////////////////
// OUTPUT FINAL BUNCH ///
/////////////////////////
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

	// Cleanup our pointers on the stack for completeness
	delete myMaterialDatabase;
	delete myBunch;
	delete myTwiss;
	delete myAccModel;
	delete myMADinterface;

    return 0;
}
