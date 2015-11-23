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
    int npart = 1E3;                          // number of particles to track
    int nturns = 200;                           // number of turns to track
 
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
	
	string input_dir = "/Thesis/data/HL_LossMap/";
	
	string output_dir = "/Build/Thesis/outputs/HL_LossMap/";
	
	bool batch = 1;
	if(batch){
		case_dir = "Merlin_1E3/";
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
	bool output_initial_bunch 	= 1;
	bool output_final_bunch 	= 0;
		if (output_initial_bunch || output_final_bunch){
			bunch_dir = (full_output_dir+"Bunch_Distn/");
			mkdir(bunch_dir.c_str(), S_IRWXU);
		}		
	bool output_fluka_database 	= 1;

///////////////////////////////////
// ACCELERATORMODEL CONSTRUCTION //
///////////////////////////////////
	cout << "MADInterface" << endl;
	
    MADInterface* myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_RF.tfs", beam_energy );
    //~ MADInterface* myMADinterface = new MADInterface( directory+input_dir+"HL_v1.2.1_C+S_original.tfs", beam_energy );

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

	string tcp_element = "TCP.C6L7.B1";    // HORIZONTAL COLLIMATOR (x)
    int tcp_element_number = myAccModel->FindElementLatticePosition(tcp_element.c_str()); 
	int start_element_number = tcp_element_number;

    cout << "Found start element TCP.C6L7 at element number " << start_element_number << endl;

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

	
	ostringstream twiss_output_file; 
    twiss_output_file << (directory+output_dir+"LatticeFunctions.dat");
    ofstream twiss_output(twiss_output_file.str().c_str());
	myTwiss->PrintTable(twiss_output);
	
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
	collimator_db->SelectImpactFactor(tcp_element, 1.0e-6);

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

    Dispersion* myDispersion = new Dispersion(myAccModel, beam_energy);
    myDispersion->FindDispersion(start_element_number);

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
    mybeam.emit_x = impact * impact * emittance * meter;
    impact =1;
    mybeam.emit_y = impact * impact * emittance * meter;
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
    ParticleBunchConstructor* myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, horizontalHaloDistribution2);

    myBunch = myBunchCtor->ConstructParticleBunch<ProtonBunch>();
    delete myBunchCtor;

    myBunch->SetMacroParticleCharge(mybeam.charge);
    
   if(output_initial_bunch){
		ostringstream bunch_output_file;
		bunch_output_file << (bunch_dir + "initial_bunch.txt");

		ofstream* bunch_output = new ofstream(bunch_output_file.str().c_str());
		myBunch->Output(*bunch_output);
		delete bunch_output;
	}

/////////////////////
// ParticleTracker //
/////////////////////

	AcceleratorModel::RingIterator beamline = myAccModel->GetRing(start_element_number);
    ParticleTracker* myParticleTracker = new ParticleTracker(beamline, myBunch);

/////////////////////////
// Collimation Process //
/////////////////////////
    CollimateProtonProcess* myCollimateProcess =new CollimateProtonProcess(2, 4);
    
    // As well as the standard output we will create a special LossMapDustbin, this will automatically sort and collate all losses
	// into a single file that we can use to plot
	LossMapDustbin* myDustbin = new LossMapDustbin;
	myCollimateProcess->SetDustbin(myDustbin);       

    // If the ScatterAtCollimator flag is true, collimation involves a full scattering simulation, if it is false, any particle to hit a collimator jaw is lost
    myCollimateProcess->ScatterAtCollimator(true);
   
    // We must assign a ScatteringModel to the CollimateProtonProcess, this will allow us to choose the type of scattering that will be used
    ScatteringModel* myScatter = new ScatteringModel;

    // MERLIN contains various ScatteringProcesses; namely the following
    // Rutherford, Elastic pn, Elastic pN, Single Diffractive, and Inelastic
    // Each of these have a MERLIN and SixTrack like version (except the inelastic)
    // As well as this there are 2 types of ionisation; simple (SixTrack like), and advanced
    // There exist 5 pre-defined set ups which may be used:
    // 0: ST,    1: ST + Adv. Ionisation,    2: ST + Adv. Elastic,    3: ST + Adv. SD,     4: MERLIN
    // Where ST = SixTrack like, Adv. = Advanced, SD = Single Diffractive, and MERLIN includes all advanced scattering

    // Below we select the simplest case (0); SixTrack like scattering and ionisation
    bool use_sixtrack_like_scattering = 0;
    if(use_sixtrack_like_scattering){
        myScatter->SetScatterType(0);
    }
    else{
        myScatter->SetScatterType(4);
    }

    // We must assign the ScatteringModel to the CollimateProtonProcess
    myCollimateProcess->SetScatteringModel(myScatter);

    // SetLossThreshold will stop the simulation if the percentage given (in this case 200%) of particles are lost
    // This case is obviously not possible and is used to make sure that the simulation is not stopped in the collimation process
    myCollimateProcess->SetLossThreshold(200.0);

    // The bin size is used to define the maximum step size in the collimator, as well as the output bin size, here it is set to 0.1m
    myCollimateProcess->SetOutputBinSize(0.1);

    // Finally we attach the process to the tracker, note that the tracker runs processes before tracking
    // The collimation process operates such that any particle outside the aperture in a collimator is collimated (scattered etc)
    // A particle outside the aperture in any other element is not collimated, instead it is immediately lost
    myParticleTracker->AddProcess(myCollimateProcess);

//////////////////
// TRACKING RUN //
//////////////////
	ostringstream particle_no_file;
        particle_no_file << pn_dir<< "No.txt";
        ofstream* particle_no_output = new ofstream(particle_no_file.str().c_str());
    if(output_turn_bunch)
		(*particle_no_output) << "0\t" << myBunch->size() << endl;
    
    // Now all we have to do is create a loop for the number of turns and use the Track() function to perform tracking   
    for (int turn=1; turn<=nturns; turn++)
    {
        // This line will give us an update of how many particles have survived after each turn
        cout << "Turn " << turn <<"\tParticle number: " << myBunch->size() << endl;

        myParticleTracker->Track(myBunch);
        
		if(output_turn_bunch)
			(*particle_no_output) << turn <<"\t" << myBunch->size() << endl;
			
        // An escape clause so that we do not needlessly track when if no particles have survived
        if( myBunch->size() <= 1 ) break;
    }
    
////////////////////////////
// OUTPUT DUSTBIN LOSSES ///
////////////////////////////
	ostringstream dustbin_file;
	dustbin_file << (full_output_dir+"HLdustbin_losses.txt");	
	ofstream* dustbin_output = new ofstream(dustbin_file.str().c_str());	
	if(!dustbin_output->good())    {
        std::cerr << "Could not open dustbin loss file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	
	// Calling Finalise() will sort and collate all losses, Output() dumps the data to the given file/stream
	myDustbin->Finalise(); 
	myDustbin->Output(dustbin_output); 
   
    // These lines tell us how many particles we tracked, how many survived, and how many were lost
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
