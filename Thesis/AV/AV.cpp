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
#include "AcceleratorModel/ApertureSurvey.h"

#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchTypes.h"
#include "BeamDynamics/ParticleTracking/Integrators/SymplecticIntegrators.h"
#include "BeamDynamics/TrackingSimulation.h"
#include "BeamDynamics/TrackingOutputASCII.h"
#include "BeamDynamics/TrackingOutputAV.h"

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

// Forward function declarations: the following functions will be used to find a specific element in the accelerator lattice
vector<AcceleratorComponent*> SortAcceleratorModel(AcceleratorModel* model);
int FindElementLatticePosition(string RequestedElement, AcceleratorModel* model);
bool SortComponent(const AcceleratorComponent* first, const AcceleratorComponent* last);

// Main function, this executable can be run with the arguments number_of_particles seed

//e.g. for 1000 particles and a seed of 356: ./test 1000 356

int main(int argc, char* argv[])
{
    int seed = (int)time(NULL);                 // seed for random number generators
    int npart = 1E3;                          // number of particles to track
    int nturns = 1;                           // number of turns to track
	bool DoTwiss = 1;							// run twiss and align to beam envelope etc?
	bool beam1 = 1;
 
    if (argc >=2){
        npart = atoi(argv[1]);
    }

    if (argc >=3){
        seed = atoi(argv[2]);
    }

// Initialise the random number generator with the seed

    RandomNG::init(seed);

    double beam_energy = 6500.0;

    cout << "npart=" << npart << " nturns=" << nturns << " beam energy = " << beam_energy << endl;

	string start_element;
	if(beam1)
	    start_element = "TCP.C6L7.B1";    // HORIZONTAL COLLIMATOR (x)
	else	
	    start_element = "TCP.C6R7.B2";    // HORIZONTAL COLLIMATOR (x)
	    
    //~ string tcp_element = "TCP.C6L7.B1";    // HORIZONTAL COLLIMATOR (x)
    //~ string start_element = "IP1";    // HORIZONTAL COLLIMATOR (x)

    // Define useful variables
    double beam_charge = 1.1e11;
    double normalized_emittance = 3.5e-6;
    double gamma = beam_energy/PhysicalConstants::ProtonMassMeV/PhysicalUnits::MeV;
	double beta = sqrt(1.0-(1.0/pow(gamma,2)));
	double emittance = normalized_emittance/(gamma*beta);
	bool output_fluka_database = 		1;
	//~ string directory = "/afs/cern.ch/user/h/harafiqu/public/MERLIN";	//lxplus harafiqu
	//~ string directory = "/home/haroon/git/Merlin/HR";				//iiaa1
	string directory = "/home/HR/Downloads/MERLIN";					//M11x	
	//~ string directory = "/afs/cern.ch/work/a/avalloni/private/MerlinforFluka/MERLIN";	//lxplus avalloni
	
	string input_dir = "/Thesis/data/AV/";
	
	//~ string output_dir = "/test2/UserSim/outputs/HL/";
	string output_dir = "/Build/Thesis/outputs/AV/";
	string batch_directory="beam2_test/";

	string full_output_dir = (directory+output_dir);
	mkdir(full_output_dir.c_str(), S_IRWXU);
	full_output_dir = (directory+output_dir+batch_directory);
	mkdir(full_output_dir.c_str(), S_IRWXU);
	
///////////////////////////////////
// ACCELERATORMODEL CONSTRUCTION //
///////////////////////////////////
	cout << "MADInterface" << endl;

    // The MADInterface class takes the MADX TFS table and the beam energy as arguments
    //    It reads the TFS table, using the keywords to build an AcceleratorModel
    //~ MADInterface* myMADinterface = new MADInterface( "/home/HR/Downloads/MERLIN/UserSim/data/TFSTable.tfs", beam_energy );
    //~ MADInterface* myMADinterface = new MADInterface( "/home/HR/Downloads/MERLIN/UserSim/data/HL_Input/twiss_hllhc_b1.tfs", beam_energy );
    //~ MADInterface* myMADinterface = new MADInterface( "/home/HR/Downloads/MERLIN/UserSim/data/HL_Input/twiss_hllhc_b1_thick.tfs", beam_energy );
    MADInterface* myMADinterface;
    if(beam1)
		myMADinterface = new MADInterface( directory+input_dir+"Twiss_6p5TeV.tfs", beam_energy );
    else
		myMADinterface = new MADInterface( directory+input_dir+"Twiss_6p5TeV_flat_top_beam2.tfs", beam_energy );
	cout << "MADInterface Done" << endl;

    // As we are only tracking for 200 turns we can choose to ignore the accelerating cavities
    // To do this we use the TreatTypeAsDrift() function, which takes an element type string as an argument, this can be done for any element
    //~ myMADinterface->TreatTypeAsDrift("RFCAVITY");
    //~ myMADinterface->TreatTypeAsDrift("SEXTUPOLE");
    //~ myMADinterface->TreatTypeAsDrift("OCTUPOLE");

    // If, as in the LHC case, the actual lattice apertures do not correspond to a single aperture per lattice element, we must not construct these apertures
    // Instead we will use an aperture file later to construct the correct apertures
    myMADinterface->ConstructApertures(false);

    // Now we can build an AcceleratorModel using MADInterface
    AcceleratorModel* myAccModel = myMADinterface->ConstructModel();    


///////////
// TWISS //
///////////


    int start_element_number_test = FindElementLatticePosition(start_element.c_str(), myAccModel);
    
    cout << "Found start element TCP.C6L7 at element number " << start_element_number_test << endl;

    // We have an AcceleratorModel, but two further items must be addressed, firstly the collimators, then the apertures
    // The LatticeFunctionTable is used to calculate the optics of the accelerator lattice, MERLIN does not use the values given in the TFS table
    LatticeFunctionTable* myTwiss = new LatticeFunctionTable(myAccModel, beam_energy);
    myTwiss->AddFunction(1,6,3);
    myTwiss->AddFunction(2,6,3);
    myTwiss->AddFunction(3,6,3);
    myTwiss->AddFunction(4,6,3);
    myTwiss->AddFunction(6,6,3);

    // Next we find the TWISS parameters
    double bscale1 = 1e-22;
    //~ double bscale1 = 1e-3;
    
    if(DoTwiss){    
		while(true)
		{
		cout << "start while(true) to scale bend path length" << endl;
			// If we are running a lattice with no RF, the TWISS parameters will not be calculated correctly
			// This is because some are calculated from using the eigenvalues of the one turn map, which is not complete with RF (i.e. longitudinal motion) switched off
			// In order to compensate for this we use the ScaleBendPath function which calls a RingDeltaT process and attaches it to the TWISS tracker
			// RingDeltaT process adjusts the ct and dp values such that the TWISS may be calculated and there are no convergence errors
			myTwiss->ScaleBendPathLength(bscale1);
			myTwiss->Calculate();

			// If Beta_x is a number (as opposed to -nan) then we have calculated the correct TWISS parameters, otherwise the loop keeps running
			if(!std::isnan(myTwiss->Value(1,1,1,0))) {break;}
			bscale1 *= 2;
			cout << "\n\ttrying bscale = " << bscale1<< endl;
		}
	}
	
	ostringstream twiss_output_file; 
    twiss_output_file << (directory+output_dir+"LatticeFunctions.dat");
    ofstream twiss_output(twiss_output_file.str().c_str());
	myTwiss->PrintTable(twiss_output);
	
///////////////////////
// Collimator set up //
///////////////////////
	cout << "Collimator Setup" << endl;   
    // We must create a new MaterialDatabase, this is essentially a standard dictionary of materials and their properties
    MaterialDatabase* myMaterialDatabase = new MaterialDatabase();
    // The collimator file is read by the CollimatorDatabase, which proceeds to set up all collimators
    CollimatorDatabase* collimator_db;
    if(beam1)
		collimator_db = new CollimatorDatabase( directory+input_dir+"Collimator_6p5TeV.txt", myMaterialDatabase,  true);
    else
		collimator_db = new CollimatorDatabase( directory+input_dir+"Collimator_6p5TeV_flat_top_beam2.txt", myMaterialDatabase,  true);
   
    // Flag to set automatic matching between beam envelope and collimator taper
    collimator_db->MatchBeamEnvelope(true);
    // Flag to set jaw alignment errors
    collimator_db->EnableJawAlignmentErrors(false);

    // If you want to set collimator jaw position or angle errors that should be done as follows
    collimator_db->SetJawPositionError(0.0 * nanometer);
    collimator_db->SetJawAngleError(0.0 * microradian);

    // Now we set the impact parameter, the beam will hit the primary collimator (TCP.C6L7.B1) at injection (in the simulation rather than the machine)
    // We use 1 micron as the standard impact parameter
    collimator_db->SelectImpactFactor(start_element, 1.0e-6);

    double impact = 6;
    // Finally we set up the collimator jaws to appropriate sizes
    try{
		if(DoTwiss)
        impact = collimator_db->ConfigureCollimators(myAccModel, emittance, emittance, myTwiss);
		else
        collimator_db->ConfigureCollimators(myAccModel);
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
    
    
     if(output_fluka_database && seed == 1){
		ostringstream fd_output_file;
		fd_output_file << (full_output_dir+"fluka_database.txt");

		ofstream* fd_output = new ofstream(fd_output_file.str().c_str());
		collimator_db->OutputFlukaDatabase(fd_output);
		delete fd_output;
	}
    delete collimator_db;
    
    //CHECK FOR COLLIMATOR APERTURES	
	vector<Collimator*> TCP;
	int siz = myAccModel->ExtractTypedElements(TCP, start_element);

	cout << "\n\t Found " << TCP.size() << " Collimators when extracting" << endl;

	Aperture *ap = (TCP[0])->GetAperture();
	if(!ap){cout << "Could not get tcp ap" << endl;	abort();}
	else{cout << "TCP aperture type = " << ap->GetApertureType() << endl;}

	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
	if(!CollimatorJaw){cout << "Could not cast" << endl;	abort();}

////////////////////////////
// Aperture Configuration //
////////////////////////////

    // The ApertureConfiguration class reads the aperture input file and assigns the appropriate apertures to each element in the lattice
    // Note that apertures do not necessarily correspond 1:1 with each element
    // If multiple apertures are defined within the scope (i.e. length) of a single element, an interpolated aperture is assigned to said element
    // The second argument to this function is the flag AllRectEllipse - if the aperture file is missing the AP_TYPE, set to true
    //~ ApertureConfiguration* myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_6p5TeV.tfs",1);     
    ApertureConfiguration* myApertureConfiguration;
    if(beam1) 
		myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_6p5TeV.tfs",1);   
    else   
		myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_6p5TeV_beam2.tfs",1);      
    
    myApertureConfiguration->ConfigureElementApertures(myAccModel);
    delete myApertureConfiguration;

// The accelerator lattice, in the form of an AcceleratorModel, is now complete
// The AcceleratorModel consists of a vector of AcceleratorComponent objects, each with it's own aperture and geometry
// Each magnet has it's own field, and each collimator has it's own material

//ApertureSurvey

	//~ Only output this file once
	//~ if(seed == 1)
	ApertureSurvey* myApertureSurvey = new ApertureSurvey(myAccModel, full_output_dir, 0.1, 5);


///////////////////
// BEAM SETTINGS //
///////////////////

    // We need to calculate the dispersion for the BeamData object
    Dispersion* myDispersion = new Dispersion(myAccModel, beam_energy);
    int start_element_number = FindElementLatticePosition(start_element.c_str(), myAccModel);
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

    // As we are tracking protons we create a proton bunch   
    ProtonBunch* myBunch;
    int node_particles = npart;

    // horizontalHaloDistribution1 is a halo in xx' plane, zero in yy'
    // horizontalHaloDistribution2 is a halo in xx' plane, gaussian in yy'
    ParticleBunchConstructor* myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, horizontalHaloDistribution2);

    myBunch = myBunchCtor->ConstructParticleBunch<ProtonBunch>();
    delete myBunchCtor;

    myBunch->SetMacroParticleCharge(mybeam.charge);

// Our bunch is now complete and ready for tracking & collimation

/////////////////////
// ParticleTracker //
/////////////////////

    // We need to create an AcceleratorModel iterator to iterate through the AcceleratorComponents in the tracker
    // The GetRing function takes the starting element number as an argument, and will return an interator that can be used for more than one turn
    AcceleratorModel::RingIterator beamline = myAccModel->GetRing(start_element_number);
    ParticleTracker* myParticleTracker = new ParticleTracker(beamline, myBunch);
    //~ myParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
    myParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());

	//~ string tof = "Tracking_output_file.txt";
	//~ string t_o_f = full_output_dir+tof;
	//~ TrackingOutputASCII* myTrackingOutputASCII = new TrackingOutputASCII(t_o_f);
	//~ myTrackingOutputASCII->SuppressUnscattered(npart+1);
	//~ myTrackingOutputASCII->output_all = 1;
	
	//~ myParticleTracker->SetOutput(myTrackingOutputASCII);
	
	ostringstream alessias_sstream;
	alessias_sstream << full_output_dir<<"Tracking_output_file"<< npart << "_" << seed << std::string(".txt");	
	string alessias_file = alessias_sstream.str().c_str();	
	
	//~ string tof = "Tracking_output_file.txt";
	//~ string t_o_f = full_output_dir+tof;
	TrackingOutputAV* myTrackingOutputAV = new TrackingOutputAV(alessias_file);
	//~ myTrackingOutputAV->SetSRange(19000, 21000);
	myTrackingOutputAV->SetSRange(0, 27000);
	myTrackingOutputAV->SetTurn(1);
	myTrackingOutputAV->output_all = 1;
	
	myParticleTracker->SetOutput(myTrackingOutputAV);

/////////////////////////
// Collimation Process //
/////////////////////////
    // Finally we create any PhysicsProcesses and assign them to the tracker
    // In this case we only need the collimation process

    // Output stream to store collimator losses
    //~ ostringstream col_output_file; 
    //~ col_output_file << (directory+output_dir+std::string("HLossMr_")) << npart << "_" << seed << std::string(".txt");
    //~ ofstream* col_output = new ofstream(col_output_file.str().c_str());
    //~ if(!col_output->good())    {
        //~ std::cerr << "Could not open collimation loss file" << std::endl;
        //~ exit(EXIT_FAILURE);
    //~ }   
   
 
    // We declare our process, every PhysicsProcess takes a priority and a mode, as we have no other processes the priority is irrelevent
    // As we have set up no modes in the CollimateParticleProcess, the mode is also irrelevent
    // We also give the output file created earlier as an argument to CollimateProtonProcess
    //~ CollimateProtonProcess* myCollimateProcess =new CollimateProtonProcess(2, 4, col_output);
    CollimateProtonProcess* myCollimateProcess =new CollimateProtonProcess(2, 4, NULL);
    
    // As well as the standard output we will create a special LossMapDustbin, this will automatically sort and collate all losses
	// into a single file that we can use to plot
	LossMapDustbin* myLossMapDustbin = new LossMapDustbin;
	myCollimateProcess->SetDustbin(myLossMapDustbin);   
	
	
	FlukaDustbin* myFlukaDustbin = new FlukaDustbin;
	myCollimateProcess->SetDustbin(myFlukaDustbin);      

    // If the ScatterAtCollimator flag is true, collimation involves a full scattering simulation, if it is false, any particle to hit a collimator jaw is lost
    myCollimateProcess->ScatterAtCollimator(true);
   
    // We must assign a ScatteringModel to the CollimateProtonProcess, this will allow us to choose the type of scattering that will be used
    ScatteringModel* myScatter = new ScatteringModel;
    if(beam1){
		myScatter->SetScatterPlot("TCP.C6L7.B1");
		myScatter->SetJawImpact("TCP.C6L7.B1");
		myScatter->SetScatterPlot("TCP.B6L7.B1");
		myScatter->SetJawImpact("TCP.D6L7.B1");
	}
	else{
		myScatter->SetScatterPlot("TCP.C6R7.B2");
		myScatter->SetJawImpact("TCP.C6R7.B2");
		myScatter->SetScatterPlot("TCP.B6R7.B2");
		myScatter->SetJawImpact("TCP.D6R7.B2");	
	}

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

    // Now all we have to do is create a loop for the number of turns and use the Track() function to perform tracking   
    for (int turn=1; turn<=nturns; turn++)
    {
        // This line will give us an update of how many particles have survived after each turn
        cout << "Turn " << turn <<"\tParticle number: " << myBunch->size() << endl;

		//~ myTrackingOutputAV->IncrementTurn();

        myParticleTracker->Track(myBunch);

        // An escape clause so that we do not needlessly track when if no particles have survived
        if( myBunch->size() <= 1 ) break;
    }

    // Flushing the output file attempts to conclude all writing to the file, then it is closed and deleted
    //~ col_output->flush();
    //~ col_output->close();
    //~ delete col_output;

	
	/*********************************************************************
	**	Output Jaw Impact
	*********************************************************************/
	myScatter->OutputJawImpact(full_output_dir);
	myScatter->OutputScatterPlot(full_output_dir);	

	/*********************************************************************
	** OUTPUT FLUKA LOSSES 
	*********************************************************************/
	ostringstream fluka_dustbin_file;
	fluka_dustbin_file << full_output_dir<<std::string("fluka_losses_")<< npart << "_" << seed << std::string(".txt");	   
	  
	ofstream* fluka_dustbin_output = new ofstream(fluka_dustbin_file.str().c_str());	
	if(!fluka_dustbin_output->good())    {
        std::cerr << "Could not open dustbin loss file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	
	myFlukaDustbin->Finalise(); 
	myFlukaDustbin->Output(fluka_dustbin_output); 
  
  
   /*********************************************************************
	** OUTPUT LOSSMAP  
	*********************************************************************/
	ostringstream dustbin_file;
	dustbin_file << full_output_dir<<"Dustbin_losses_"<< npart << "_" << seed << std::string(".txt");	
	ofstream* dustbin_output = new ofstream(dustbin_file.str().c_str());	
	if(!dustbin_output->good())    {
        std::cerr << "Could not open dustbin loss file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	
	myLossMapDustbin->Finalise(); 
	myLossMapDustbin->Output(dustbin_output); 
   
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



// The following functions will be used to find a specific element in the accelerator lattice
// Firstly a simple comparison function to establish which element comes first in the lattice depending on its (s) position

bool SortComponent(const AcceleratorComponent* first, const AcceleratorComponent* last){
    return (first->GetComponentLatticePosition() < last->GetComponentLatticePosition());
}

// Next a function that returns a vector of AcceleratorComponents that have been sorted in order of position using the previous function

vector<AcceleratorComponent*> SortAcceleratorModel(AcceleratorModel* model){
    vector<AcceleratorComponent*> elements;

    model->ExtractTypedElements(elements,"*"); //This line extracts all elements as the wildcard * is used

    //Now sort the elements in the appropriate lattice order
    sort(elements.begin(), elements.end(),SortComponent);
    return elements;
}

// Finally a function to find the position of an element in the lattice using the previous function

int FindElementLatticePosition(string RequestedElement, AcceleratorModel* model){
    vector<AcceleratorComponent*> elements = SortAcceleratorModel(model);
    size_t nelm = elements.size();
    for(size_t n=0; n<nelm; n++){
        if(elements[n]->GetName() == RequestedElement){
            cout << "Found " << RequestedElement << " at " << n << " of " << nelm << endl;
            return n;
        }
    }
    return 0;
}
