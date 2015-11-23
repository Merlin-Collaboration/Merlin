#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <unistd.h>

#include "AcceleratorModel/Components.h"
#include "AcceleratorModel/Apertures/CollimatorAperture.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"

#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchTypes.h"

#include "Collimators/CollimateParticleProcess.h"
#include "Collimators/CollimateProtonProcess.h"
#include "Collimators/MaterialDatabase.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"

#include "Random/RandomNG.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;


int main(int argc, char* argv[])
{
	int seed = (int)time(NULL);
	if (argc >=2)
    {
        seed = atoi(argv[1]);
    }	

	cout << "Seed: " << seed << endl;
	RandomNG::init(seed);
	/*********************************************************************
	**	GENERAL SETTINGS
	*********************************************************************/
	//Loss_Map or Merged Collimation
    bool Loss_Map 				= 0;
		if(Loss_Map){std::cout << "LOSSMAP Collimation (ProtonBunch)" << std::endl;}
		else{std::cout << "MERGED Collimation (based on HEL code)" << std::endl;}
    bool output_initial_bunch 	= 0;
    bool output_final_bunch 	= 0;
	
	//Beam energy (GeV) 7000,3500,450 etc
	//double beam_energy = 7000.0;
	const double beam_energy = 7000.0;

	//Number of particles
	const int npart = 1E4;
	//~ const int npart = 100;
	
	/*********************************************************************
	**	ACCELERATOR MODEL LOADING
	*********************************************************************/

  	MaterialDatabase* mat = new MaterialDatabase();
  	//~ Material* CollimatorMaterial = mat->FindMaterial("W");
  	Material* CollimatorMaterial = mat->FindMaterial("Cu");

	AcceleratorModelConstructor* construct = new AcceleratorModelConstructor();	
	construct->NewModel();
	double length = 1;
	Collimator* TestCol = new Collimator("TestCollimator",length);
	TestCol->SetMaterial(CollimatorMaterial);

	CollimatorAperture* app=new CollimatorAperture(2, 2, 0, CollimatorMaterial, length, 0, 0);
	app->SetExitWidth(app->GetFullEntranceWidth());      //Horizontal
	app->SetExitHeight(app->GetFullEntranceHeight());    //Vertical

	//Set the aperture for collimation
	TestCol->SetAperture(app);

	construct->AppendComponent(*TestCol);
	
	AcceleratorModel* model = construct->GetModel();

	/*********************************************************************
	**      BEAM SETTINGS
	*********************************************************************/
	ProtonBunch* myBunch = new ProtonBunch(beam_energy,1);
	Particle p(0);
	p.y() = 1.0 + 1e-6;
	for(int i=0; i<npart; i++)
	{
		p.id() = i;
		myBunch->AddParticle(p);
	}
	/*********************************************************************
	**	Output Initial Bunch
	*********************************************************************/
	if(output_initial_bunch)
	{
		ostringstream bunch_output_file;
		if(Loss_Map)
			bunch_output_file << "Thesis/outputs/LM_ST_initial.txt";
		else
			bunch_output_file << "Thesis/outputs/HEL_ST_initial.txt";
		

		ofstream* bunch_output = new ofstream(bunch_output_file.str().c_str());
		myBunch->Output(*bunch_output);
		delete bunch_output;
	 }
	 
	/*********************************************************************
	**	PARTICLE TRACKER
	*********************************************************************/
	AcceleratorModel::RingIterator bline = model->GetRing();
	ParticleTracker* tracker = new ParticleTracker(bline,myBunch);

	/*********************************************************************
	**	COLLIMATION SETTINGS
	*********************************************************************/
		
	CollimateProtonProcess* myCollimateProcess = new CollimateProtonProcess(2, 4, NULL);

	LossMapDustbin* myDustbin = new LossMapDustbin(tencm);
	myCollimateProcess->SetDustbin(myDustbin);       
	myCollimateProcess->ScatterAtCollimator(true);
   
	ScatteringModel* myScatter = new ScatteringModel;

	bool use_sixtrack_like_scattering = 0;
	if(use_sixtrack_like_scattering){
		myScatter->SetScatterType(0);
	}
	else{
		myScatter->SetScatterType(4);
	}
	
	myScatter->SetScatterPlot("TestCollimator");
	myScatter->SetJawImpact("TestCollimator");

	myCollimateProcess->SetScatteringModel(myScatter);

	myCollimateProcess->SetLossThreshold(200.0);

	myCollimateProcess->SetOutputBinSize(0.1);
	//~ myCollimateProcess->SetOutputBinSize(length);

	tracker->AddProcess(myCollimateProcess);


	/*********************************************************************
	**	Tracking
	*********************************************************************/

	cout << "Tracking" << endl;
	tracker->Track(myBunch);
	cout << "Finished.\tParticle number: " << myBunch->size() << endl;
	cout << "npart: " << npart << endl;
	cout << "left: " << myBunch->size() << endl;
	cout << "absorbed: " << npart - myBunch->size() << endl;

	if (0){
		ostringstream bunch_output_file_out;
		bunch_output_file_out << "cu50_test_bunch_out.txt";
		ofstream* bunch_output_out = new ofstream(bunch_output_file_out.str().c_str());
		*bunch_output_out << "#T0 P0 x xp y yp ct dp" << endl;
		myBunch->Output(*bunch_output_out);
	}
	
	/*********************************************************************
	**	Output Scatter Plot
	*********************************************************************/
	//~ ostringstream scatter_file;
	//~ scatter_file << "Thesis/outputs/ScatterPlot.txt";
	//~ ofstream* scatter_output = new ofstream(scatter_file.str().c_str());
			
	myScatter->OutputScatterPlot("Thesis/outputs/");
	
	//~ delete scatter_output;
	
	/*********************************************************************
	**	Output Jaw Impact
	*********************************************************************/
	//~ ostringstream impact_file;
	//~ impact_file << "Thesis/outputs/JawImpact.txt";
	//~ ofstream* impact_output = new ofstream(impact_file.str().c_str());
			
	myScatter->OutputJawImpact("Thesis/outputs/");
	
	//~ delete impact_output;
	
	/*********************************************************************
	**	Output Final Bunch
	*********************************************************************/
	if(output_final_bunch){
		ostringstream bunch_output_file2;
		if(Loss_Map)
			bunch_output_file2 << "Thesis/outputs/LM_M_final.txt";
		else
			bunch_output_file2 << "Thesis/outputs/HEL_M_final.txt";

		ofstream* bunch_output2 = new ofstream(bunch_output_file2.str().c_str());
		myBunch->Output(*bunch_output2);
		delete bunch_output2;
	 }

	delete myBunch;

	return 0;
}
