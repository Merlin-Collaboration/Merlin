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
    int npart = 100;                            // number of particles to track
    int nturns = 1;                           // number of turns to track
	bool DoTwiss = 1;							// run twiss and align to beam envelope etc?
	
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

//HEL SETTINGS		
	bool nominal 				= 0;			// Use nominal or HL LHC		
	bool LHC_HEL				= 1;			// Tevatron or LHC HEL parameters
	if(!nominal){LHC_HEL = 1;}
	
	bool hel_on 				= 0; 			// Hollow electron lens process on?
		bool DCon							= 0;
		bool ACon							= 0;
		if(ACon){DCon=0;}
		bool Turnskipon						= 0;
		if(Turnskipon){ACon=0; DCon=0;}
		bool Diffusiveon					= 1;
		if(Diffusiveon){ACon=0; Turnskipon=0; DCon=0;}

	bool output_fluka_database	= 1;
	
//OUTPUTS   
    //Working Directory
	//~ string directory = "/afs/cern.ch/user/h/harafiqu/public/MERLIN";	//lxplus harafiqu
	//~ string directory = "/home/haroon/git/Merlin";						//iiaa1
	string directory = "/home/HR/Downloads/MERLIN";							//M11x	
	
	string input_dir, case_dir, pn_dir, bunch_dir;
		
	input_dir = "/Thesis/data/JawImpact/";
	
	string output_dir = "/Build/Thesis/outputs/JawImpact/";
	
	bool batch = 1;
	if(batch){
		// no HEL
		//~ case_dir = "Nom_1E6/";
		//~ case_dir = "HL_1E6/";		
		// Diffusive HEL
		//~ case_dir = "Nom_LHC_1E7_Diff/";
		//~ case_dir = "Nom_Tev_1E7_Diff/";
		//~ case_dir = "HL_LHC_1E7_Diff/";
		//~ case_dir = "Symplectic/";
		case_dir = "Transport_10/";
	}
	
	string full_output_dir = (directory+output_dir);
	mkdir(full_output_dir.c_str(), S_IRWXU);
	
	if(batch){
		full_output_dir = (directory+output_dir+case_dir);
		mkdir(full_output_dir.c_str(), S_IRWXU);
	}
	
	bool output_initial_bunch 	= 1;
	bool output_final_bunch 	= 0;
		if (output_initial_bunch || output_final_bunch){
			bunch_dir = (full_output_dir+"Bunch_Distn/");
			mkdir(bunch_dir.c_str(), S_IRWXU);
		}
	
	bool output_turn_bunch		= 0;
		if(output_turn_bunch){
			pn_dir = (full_output_dir+"ParticleNo/");
			mkdir(pn_dir.c_str(), S_IRWXU);
		}
	bool output_scatter_plot 	= 1;
	bool output_jaw_impact		= 1;
	
///////////////////////////////////
// ACCELERATORMODEL CONSTRUCTION //
///////////////////////////////////
	cout << "MADInterface" << endl;
	MADInterface* myMADinterface;
	
	if(nominal) //old modified 7TeV+HEL lattice
	myMADinterface = new MADInterface( directory+input_dir+"Nominal+HEL.tfs", beam_energy );
	else if(!nominal)	//new HL+HEL lattice
	myMADinterface = new MADInterface( directory+input_dir+"HLv1.2.0+HEL.tfs", beam_energy );

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
    
	string ip1_element = "IP1";    // IP1   
    int ip1_element_number = myAccModel->FindElementLatticePosition(ip1_element.c_str()); 
    
    //~ int start_element_number = tcp_element_number;
    int start_element_number = ip1_element_number;
    
    //~ cout << "Found start element TCP.C6L7 at element number " << start_element_number << endl;
    cout << "Found start element at element number " << start_element_number << endl;

    LatticeFunctionTable* myTwiss = new LatticeFunctionTable(myAccModel, beam_energy);
    myTwiss->AddFunction(1,6,3);
    myTwiss->AddFunction(2,6,3);
    myTwiss->AddFunction(3,6,3);
    myTwiss->AddFunction(4,6,3);
    myTwiss->AddFunction(6,6,3);

    double bscale1 = 1e-22;
    
    if(DoTwiss){    
		while(true)
		{
		cout << "start while(true) to scale bend path length" << endl;
			myTwiss->ScaleBendPathLength(bscale1);
			myTwiss->Calculate();

			if(!std::isnan(myTwiss->Value(1,1,1,0))) {break;}
			bscale1 *= 2;
		}
	}
	
	ostringstream twiss_output_file; 
    twiss_output_file << (directory+output_dir+"LatticeFunctions.dat");
    ofstream twiss_output(twiss_output_file.str().c_str());
	myTwiss->PrintTable(twiss_output);

////////////////////////////
// Aperture Configuration //
////////////////////////////
	ApertureConfiguration* myApertureConfiguration;
    
    if(nominal) //old modified 7TeV+HEL lattice      
    myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_Nominal+HEL.tfs",0); 
    else if(!nominal)	//new HL+HEL lattice     
    myApertureConfiguration = new ApertureConfiguration(directory+input_dir+"Aperture_HLv1.2.0.tfs",1);      
    
    myApertureConfiguration->ConfigureElementApertures(myAccModel);
    delete myApertureConfiguration;
    	
///////////////////////
// Collimator set up //
///////////////////////
	cout << "Collimator Setup" << endl;   
    MaterialDatabase* myMaterialDatabase = new MaterialDatabase();
    CollimatorDatabase* collimator_db;    
    
    //Need to use proper collimator settings
    if(nominal) //old modified 7TeV+HEL lattice
	collimator_db = new CollimatorDatabase( directory+input_dir+"ColDB_Nominal.txt", myMaterialDatabase,  true);    
    else if(!nominal)	//new HL+HEL lattice
	collimator_db = new CollimatorDatabase( directory+input_dir+"ColDB_HLv1.2.0.txt", myMaterialDatabase,  true);
   
   
    collimator_db->MatchBeamEnvelope(true);
    collimator_db->EnableJawAlignmentErrors(false);

    collimator_db->SetJawPositionError(0.0 * nanometer);
    collimator_db->SetJawAngleError(0.0 * microradian);

    collimator_db->SelectImpactFactor(tcp_element, 1.0e-6);

    double impact = 6;
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
		if(!ap)
			{cout << "Could not get tcp ap" << endl;	abort();}
		else{cout << "TCP aperture type = " << ap->GetApertureType() << endl;}

		 CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
		if(!CollimatorJaw)
			{cout << "Could not cast" << endl;	abort();}
			
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

    // As we are tracking protons we create a proton bunch   
    ProtonBunch* myBunch;
    int node_particles = npart;

    // horizontalHaloDistribution1 is a halo in xx' plane, zero in yy'
    // horizontalHaloDistribution2 is a halo in xx' plane, gaussian in yy'
    //~ ParticleBunchConstructor* myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, horizontalHaloDistribution2);
    ParticleBunchConstructor* myBunchCtor = new ParticleBunchConstructor(mybeam, node_particles, LHCDistn);

    myBunch = myBunchCtor->ConstructParticleBunch<ProtonBunch>();
    delete myBunchCtor;

    myBunch->SetMacroParticleCharge(mybeam.charge);

// Our bunch is now complete and ready for tracking & collimation
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

    AcceleratorModel::Beamline prebeamline = myAccModel->GetBeamline(ip1_element_number, tcp_element_number-1);
    AcceleratorModel::RingIterator beamline = myAccModel->GetRing(tcp_element_number);
    
    ParticleTracker* myParticleTracker = new ParticleTracker(beamline, myBunch);
	ParticleTracker* myPreParticleTracker = new ParticleTracker(prebeamline, myBunch);
        
	myPreParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());	
	myParticleTracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());
	
	//~ myPreParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
	//~ myParticleTracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
		

/////////////////////////
// Collimation Process //
/////////////////////////

    CollimateProtonProcess* myCollimateProcess =new CollimateProtonProcess(2, 4, NULL);
    
 	LossMapDustbin* myLossMapDustbin = new LossMapDustbin;	 
		myCollimateProcess->SetDustbin(myLossMapDustbin); 
		
	FlukaDustbin* myFlukaDustbin = new FlukaDustbin;	
		myCollimateProcess->SetDustbin(myFlukaDustbin);  
	 
    myCollimateProcess->ScatterAtCollimator(true);
   
    ScatteringModel* myScatter = new ScatteringModel;

    bool use_sixtrack_like_scattering = 0;
    if(use_sixtrack_like_scattering){
        myScatter->SetScatterType(0);
    }
    else{
        myScatter->SetScatterType(4);
    }    
    
	myScatter->SetScatterPlot("TCP.C6L7.B1");
	myScatter->SetJawImpact("TCP.C6L7.B1");
	//~ myScatter->SetScatterPlot("Collimator.TCP.C6L7.B1");
	//~ myScatter->SetJawImpact("Collimator.TCP.C6L7.B1");

    myCollimateProcess->SetScatteringModel(myScatter);

    myCollimateProcess->SetLossThreshold(200.0);

    myCollimateProcess->SetOutputBinSize(0.1);

    myParticleTracker->AddProcess(myCollimateProcess);
    myPreParticleTracker->AddProcess(myCollimateProcess);
    
/////////////////
// HEL Process //
/////////////////	
	if(hel_on){
		
		// HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e);
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
			myHELProcess->SetTurnskip(31);
			myHELProcess->SetOpMode(Turnskip);
		}	
		myParticleTracker->AddProcess(myHELProcess);	
	}
//////////////////
// TRACKING RUN //
//////////////////
	ostringstream particle_no_file;
        particle_no_file << pn_dir<< "No.txt";
        ofstream* particle_no_output = new ofstream(particle_no_file.str().c_str());

	//PRETRACK
	myPreParticleTracker->Track(myBunch);
        
	if(output_turn_bunch)
		(*particle_no_output) << "0\t" << myBunch->size() << endl;
		
    for (int turn=1; turn<=nturns; turn++)
    {
        cout << "Turn " << turn <<"\tParticle number: " << myBunch->size() << endl;

        myParticleTracker->Track(myBunch);
		
		if(output_turn_bunch)
			(*particle_no_output) << turn <<"\t" << myBunch->size() << endl;				
	
	    if( myBunch->size() <= 1 ) break;
    }
    
////////////////////////
//  END TRACKING RUN  //
////////////////////////
	if(output_scatter_plot){    
		/*********************************************************************
		**	Output Scatter Plot
		*********************************************************************/
		myScatter->OutputScatterPlot(full_output_dir);		
	}
	//~ if(output_jaw_impact && myScatter->StoredJawImpactData.size() != 0){
	if(output_jaw_impact){
		/*********************************************************************
		**	Output Jaw Impact
		********************************************************************/
		myScatter->OutputJawImpact(full_output_dir);
	}
	/*********************************************************************
	** OUTPUT FLUKA DUSTBIN 
	*********************************************************************/
	ostringstream fluka_dustbin_file;
	fluka_dustbin_file << (full_output_dir+"fluka_losses.txt");	
	ofstream* fluka_dustbin_output = new ofstream(fluka_dustbin_file.str().c_str());	
	if(!fluka_dustbin_output->good())    {
        std::cerr << "Could not open dustbin loss file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	
	myFlukaDustbin->Finalise(); 
	myFlukaDustbin->Output(fluka_dustbin_output); 
	
	/*********************************************************************
	** OUTPUT LOSSMAP DUSTBIN 
	*********************************************************************/
	ostringstream dustbin_file;
	dustbin_file << (full_output_dir+"dustbin_losses.txt");	
	ofstream* dustbin_output = new ofstream(dustbin_file.str().c_str());	
	if(!dustbin_output->good())    {
        std::cerr << "Could not open dustbin loss file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	
	myLossMapDustbin->Finalise(); 
	myLossMapDustbin->Output(dustbin_output); 
   
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

