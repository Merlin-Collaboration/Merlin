#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <map>
#include <set>
#include <ctime>


#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchTypes.h"
#include "MADInterface/MADInterface.h"
#include "Random/RandomNG.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
//~ #include "Collimators/CollimateParticleProcess.h"
#include "Collimators/CollimateProtonProcess.h"
#include "Collimators/ScatteringProcess.h"
#include "Collimators/ScatteringModel.h"
#include "Collimators/CollimatorDatabase.h"
#include "Collimators/MaterialDatabase.h"
#include "Collimators/ApertureConfiguration.h"
#include "RingDynamics/Dispersion.h"

using namespace std;
using namespace PhysicalUnits;

/*
 * Note this test should fail occasionally. Actually it fails quite often.
 * By default it only uses 1000 particles, which is too few. 10k or 100k
 * are really needed for a more reliable test, but this makes it slow.
 *
 * Compute the loss map for the nominal LHC lattice, and compare with
 * a pre-computed version.
 * 
 * arguments:
 *    collimation_test nparticles seed
 *
 * The number of bins deviating by more than N sigma are recorded
 * and if the are significantly more that expected (by a normal
 * distribution), the test fails.
 *
 * If the test fails run it a few times. If it consistently fails then
 * there is an actual issue. Bins out by more than 3 sigma are
 * displayed, repeated failure of the same bin would be suspicious.
 *
 * The time (seconds since epoch) is used as the seed unless an int is 
 * passed as the first argument.
 * 
 *
 */


enum loss_map_mode_t {HORIZONTAL_LOSS, VERTICAL_LOSS};

bool SortComponent(const AcceleratorComponent* first, const AcceleratorComponent* last)
{
	return (first->GetComponentLatticePosition() < last->GetComponentLatticePosition());
}


vector<AcceleratorComponent*> SortAcceleratorModel(AcceleratorModel* model)
{
	vector<AcceleratorComponent*> elements;
	model->ExtractTypedElements(elements,"*");

    //Now sort the elements in the appropriate lattice order
	sort(elements.begin(), elements.end(),SortComponent);
	return elements;
}

int FindElementLatticePosition(string RequestedElement, AcceleratorModel* model)
{
	vector<AcceleratorComponent*> elements = SortAcceleratorModel(model);
	size_t nelm = elements.size();
	for(size_t n=0; n<nelm; n++)
	{
		if(elements[n]->GetName() == RequestedElement)
		{
			cout << "Found " << RequestedElement << " at " << n << " of " << nelm << endl;
			return n;
		}
	}
	return 0; 
}

int main(int argc, char* argv[])
{
	int seed = (int)time(NULL);
	int npart = 1000;
	int nturns = 200;
	
	//Loss_Map or Merged Collimation
    bool Loss_Map 				= 0;

	if (argc >=2)
	{
		npart = atoi(argv[1]);
	}

	if (argc >=3)
	{
		seed = atoi(argv[2]);
	}

	cout << "Random Seed: " << seed << endl;
	RandomNG::init(seed);

	// find directories 
	string log_dir = "";
	string input_data_dir = "";
	string result_dir = "";

	string paths[] = {"../", "", "MerlinTests/"};
	for (size_t i=0; i<3; i++){
		ifstream test_file;
		test_file.open((paths[i]+"data/collimator.7.0.sigma").c_str());
		if (test_file){
			test_file.close();
			input_data_dir = paths[i]+"data/";
			log_dir = paths[i]+"outputs/";
			result_dir = paths[i]+"outputs/";
			break;
		}
	}
	if (log_dir == ""){
		cout << "Could not find data directory. Try running from cmake build directory, or the executable directory" << endl;
		exit(1);
	}

	double beam_energy = 7000.0;

	cout << "npart=" << npart << " nturns=" << nturns << endl;
	
	loss_map_mode_t loss_map_mode = HORIZONTAL_LOSS;
//	loss_map_mode_t loss_map_mode = VERTICAL_LOSS;

	string start_element;
	switch (loss_map_mode){
		case HORIZONTAL_LOSS:
			start_element = "TCP.C6L7.B1";	//HORIZONTAL COLLIMATOR (x)
			break;
		case VERTICAL_LOSS:
			start_element = "TCP.D6L7.B1";	//VERTICAL COLLIMATOR (y)
			break;
	}

	double beamcharge = 1.1e11;
	double normalized_emittance = 3.5e-6;
	double gamma = beam_energy/PhysicalConstants::ProtonMassMeV/PhysicalUnits::MeV;
	double beta = sqrt(1.0-(1.0/pow(gamma,2)));
	double emittance = normalized_emittance/(gamma*beta);

  	MaterialDatabase* mat = new MaterialDatabase();
	CollimatorDatabase* collimator_db = new CollimatorDatabase(input_data_dir+"collimator.7.0.sigma", mat, true);
	

	//	ACCELERATOR MODEL LOADING
	MADInterface* myMADinterface;
	myMADinterface = new MADInterface(input_data_dir+"twiss.7.0tev.b1_new.tfs",beam_energy);
	myMADinterface->TreatTypeAsDrift("RFCAVITY");
	myMADinterface->ConstructApertures(false);

	//Build accelerator model
	AcceleratorModel* model = myMADinterface->ConstructModel();

    LatticeFunctionTable* twiss = new LatticeFunctionTable(model,beam_energy);
	twiss->AddFunction(1,6,3);
	twiss->AddFunction(2,6,3);
	twiss->AddFunction(3,6,3);
	twiss->AddFunction(4,6,3);
	twiss->AddFunction(6,6,3);

	// find twiss
	double bscale1 = 1e-22;
	while(true)
	{
		twiss->ScaleBendPathLength(bscale1);
		twiss->Calculate();
		if(!std::isnan(twiss->Value(1,1,1,0)))
		{
			break;
		}
		bscale1 *= 2;
	}
	
	// FLAG to set an automatically matching between beam envelope and collimator taper
	collimator_db->MatchBeamEnvelope(false);
	collimator_db->EnableJawAlignmentErrors(false);

	//Collimator rms error on gap size 0.1 sigma, jaw angle error with respect the beam envelope 200 microradiant
	collimator_db->SetJawPositionError(0.0 * nanometer);// This is actually the variance of the error alignment of 0.1 sigma i.e. 35.3 micronmeter.... it should be in nm^2
	collimator_db->SetJawAngleError(0.0 * microradian);// rms error 200 microradiant
	collimator_db->SelectImpactFactor(start_element, 1.0e-6);
	
	double impact;
    //Setup the collimator jaws to appropriate sizes and
	try
	{
		impact = collimator_db->ConfigureCollimators(model, emittance, emittance,twiss);
	}
	catch(exception& e)
	{
		std::cout << "Exception caught: " << e.what() << std::endl;
		exit(1);
	}
	if(std::isnan(impact))
	{
		cerr << "Impact is nan" << endl;
		exit(1);
	}
    cout << "Impact factor number of sigmas: " << impact << endl;
	delete collimator_db;

	ApertureConfiguration* apc = new ApertureConfiguration(input_data_dir+"LHCB1Aperture.tfs");
	
    apc->ConfigureElementApertures(model);
	cout << "aperture load finished" << endl<< "start twiss" << endl;;
	delete apc;

	//Calculate Dispersion
	Dispersion* disp = new Dispersion(model,beam_energy);
	int start_element_number = FindElementLatticePosition(start_element.c_str(), model);
	disp->FindDispersion(start_element_number);

    //      BEAM SETTINGS
	BeamData mybeam;

	//Default values are 0.0
	//The charge of the particles in the beam.
	//  <0 for electrons, >0 for positrons/protons.
	mybeam.charge = beamcharge/npart;
	mybeam.p0 = beam_energy;
	mybeam.beta_x = twiss->Value(1,1,1,start_element_number)*meter;
	mybeam.beta_y = twiss->Value(3,3,2,start_element_number)*meter;
	mybeam.alpha_x = -twiss->Value(1,2,1,start_element_number);
	mybeam.alpha_y = -twiss->Value(3,4,2,start_element_number);

	//Dispersion
	mybeam.Dx=disp->Dx;
	mybeam.Dy=disp->Dy;
	//mybeam.Dy=0;
	mybeam.Dxp=disp->Dxp;
	mybeam.Dyp=disp->Dyp;

	mybeam.emit_x = impact * impact * emittance * meter;
	impact=1;
	mybeam.emit_y = impact * impact * emittance * meter;
	mybeam.sig_z = 0.0;

	//Beam centroid
	mybeam.x0=twiss->Value(1,0,0,start_element_number);
	mybeam.xp0=twiss->Value(2,0,0,start_element_number);
	mybeam.y0=twiss->Value(3,0,0,start_element_number);
	mybeam.yp0=twiss->Value(4,0,0,start_element_number);
	mybeam.ct0=twiss->Value(5,0,0,start_element_number);

	mybeam.sig_dp = 0.0;

	//X-Y coupling
	mybeam.c_xy=0.0;
	mybeam.c_xyp=0.0;
	mybeam.c_xpy=0.0;
	mybeam.c_xpyp=0.0;

	delete disp;
	
	//	BUNCH SETTINGS
	ProtonBunch* myBunch;
	int node_particles = npart;

	ParticleBunchConstructor* constructor;
	// dist 1 is halo in 1 plane, zero in other
	// dist 2 is halo in 1 plane, gauss in other
	switch (loss_map_mode){
		case HORIZONTAL_LOSS:
			constructor = new ParticleBunchConstructor(mybeam,node_particles,horizontalHaloDistribution2);
			break;
		case VERTICAL_LOSS:
			constructor = new ParticleBunchConstructor(mybeam,node_particles,verticalHaloDistribution2);
			break;
	}

	myBunch = constructor->ConstructParticleBunch<ProtonBunch>();
	delete constructor;

	myBunch->SetMacroParticleCharge(mybeam.charge);
	
	if(Loss_Map){
		//~ myBunch->EnableScatteringPhysics(ProtonBunch::Merlin);
		myBunch->EnableScatteringPhysics(ProtonBunch::SixTrack);
	}
	
	//	PARTICLE TRACKER
	AcceleratorModel::RingIterator bline = model->GetRing(start_element_number);
	ParticleTracker* tracker = new ParticleTracker(bline,myBunch);
	tracker->SetLogStream(std::cout);
	
	//	COLLIMATION SETTINGS	
	
	//Output stream for collimator losses
	ostringstream col_output_file;
	col_output_file << result_dir << "/Loss.txt";
	
	if (1){ // disable tracking
		ofstream* col_output = new ofstream(col_output_file.str().c_str());
	if(!col_output->good())
	{
		std::cerr << "Could not open collimation loss file" << std::endl;
		exit(EXIT_FAILURE);
	}	
	
	if(Loss_Map){		
		CollimateParticleProcess* myCollimateProcess;		
		
		myCollimateProcess=new CollimateParticleProcess(2,4,col_output);		
		myCollimateProcess->ScatterAtCollimator(true);

		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);
		tracker->AddProcess(myCollimateProcess);				
	}
	else{
		CollimateProtonProcess* myCollimateProcess;		
		
		myCollimateProcess=new CollimateProtonProcess(2,4,col_output);		
		myCollimateProcess->ScatterAtCollimator(true);
		
		ScatteringModel* myScatter = new ScatteringModel;
		//~ myScatter->SetScatterType(4);
		myScatter->SetScatterType(0);

		myCollimateProcess->SetScatteringModel(myScatter);

		myCollimateProcess->SetLossThreshold(200.0);
		myCollimateProcess->SetOutputBinSize(0.1);
		tracker->AddProcess(myCollimateProcess);
		
	}	
		 
	//	TRACKING RUN
	for (int turn=1; turn<=nturns; turn++)
	{
		cout << "Turn " << turn <<"\tParticle number: " << myBunch->size() << endl;
		tracker->Track(myBunch);
		if(myBunch->size() <= 1) break;

	}
	col_output->flush();
	col_output->close();
	delete col_output;
	}
	
	cout << "npart: " << npart << endl;
	cout << "left: " << myBunch->size() << endl;
	cout << "absorbed: " << npart - myBunch->size() << endl;

	delete mat;
	delete myBunch;
	delete twiss;
	delete model;
	delete myMADinterface;

	// read losses back in
	ifstream* col_output2 = new ifstream(col_output_file.str().c_str());

	set<string> all_loss_elems;
	map<string, int> loss_map;
	string el_name, el_pos, el_spos, hit_id;
	double el_len, el_hits;

	while (col_output2->good()){
		*col_output2 >> el_name >> el_pos >> el_len >> el_hits >> el_spos;
		hit_id = el_name + "_" + el_pos;
		//if(!loss_map.count(hit_id)){ loss_map[hit_id] = 0;}
		//cout << "~" << hit_id << " " <<  el_hits << endl;
		loss_map[hit_id] += el_hits;
		all_loss_elems.insert(hit_id);
	}

	// load test data
	ostringstream nom_input_file;
	nom_input_file << input_data_dir << "/lhc_collimation_test_7tev_nom.dat";

	el_hits = 7;
	map<string, double> loss_map_nom;
	ifstream* nom_lossmap_in = new ifstream(nom_input_file.str().c_str());
	string line;
	while (getline(*nom_lossmap_in, line)){
		if (line == "") continue;
		if (line[0] == '#')continue;
		istringstream liness(line);
		liness >> hit_id >> el_hits;
		//cout << "@" << hit_id << " " <<el_hits*npart << endl;;
		//if(!loss_map_nom.count(hit_id)){ loss_map_nom[hit_id] = 0;}
		loss_map_nom[hit_id] += el_hits * npart;
		all_loss_elems.insert(hit_id);
	}
	
	// Calculate sigma from sqrt(<n>)
	// Count number of bins out by more than n-sigma
	const int ns = 4;
	// https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule, could use erf() in c++11
	const double tails[] = {0, 0.682689492137086, 0.954499736103642, 0.997300203936740,
	                           0.999936657516334, 0.999999426696856, 0.999999998026825};
	int over_ns[ns+1] = {0};
	int nbins = 0;

	cout << "**" << endl;
	double sum_x1 = 0, sum_x2 = 0;
	for(set<string>::iterator it=all_loss_elems.begin(); it!=all_loss_elems.end(); ++it){
		double x1 = loss_map[*it];
		double x2 = loss_map_nom[*it];
		sum_x1 += x1;
		sum_x2 += x2;

		// small bins throw the result, so ignore them until proper poisson stats implemented FIXME
		if (x2 < 5) {continue;}
		nbins++;

		//double exp_err_x = max(sqrt(x2), 0.5);
		double exp_err_x = sqrt(x2);
        double diff_x = x1 - x2;

		for(int j=1; j<ns+1; j++){
			if (fabs(diff_x)/exp_err_x > j) {
				over_ns[j]++;
				if (j>=2) {
					cout << "big error in x # bin, merlin, expected, error/sigma" << endl;
					cout << *it << " " << x1 << " " << x2 << "+-" << exp_err_x << " " << diff_x/exp_err_x << endl;
				}
			}	
		}
		//cout << *it << " " << loss_map[*it] << " " << loss_map_nom[*it] << endl;
	}
	cout << "Number of particles used "<< npart << endl;
	
	cout << "Number of bins with sufficient expected losses "<<nbins << endl;

	if (npart < 10000){
		cout << "Warning low statistics, it is recommended to use more than 10000 particles" << endl;
	}
	
	// put limit on fraction of bins out by more than n-sigma
	// bad if 30% more than expected
	
	bool bad = 0;
	double tol = 1.3;
	cout << endl << "Bins out by more than";
	for(int j=1; j<ns+1; j++){
		cout << endl << "    "<< j <<" sigma: "<< over_ns[j] << " ("<< (double)over_ns[j]/(nbins+1)*100  << " %)";
		if ((double)over_ns[j]/(nbins+1) > (1-tails[j])*tol) {cout << "** bad **"; bad=1;}
	}
	cout << endl << endl;

	cout << "Total losses " << sum_x1 << "  expected " << sum_x2 << endl;
	double sum_dev = fabs(sum_x1-sum_x2) / sqrt(sum_x2);
	cout << sum_dev <<  " sig deviation";
	if (sum_dev > 2) {bad += 1; cout << " ** bad **";}
	cout << endl;
	

	return bad;
}
