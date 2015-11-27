#ifndef ScatteringModel_h
#define ScatteringModel_h 1

#include <iostream>
#include <cmath>
#include <map>

#include "merlin_config.h"

#include "BeamModel/PSvector.h"

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

#include "Collimators/Material.h"
#include "Collimators/ScatteringProcess.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

namespace Collimation {
	
struct JawImpactData{
	int turn;
	int ID;
	double x;
	double xp;
	double y;
	double yp;
	double ct;
	double dp;	
	string name;
};

struct ScatterPlotData{
	int turn;
	int ID;
	double x;
	double xp;
	double y;
	double yp;
	double z;
	string name;
	//~ double ct;
	//~ double dp;		
	
	inline bool operator==(const ScatterPlotData& rhs){
		if( (this->z != rhs.z) )
		{ return 0;}
		else {return 1;}
	}
	
	inline bool operator>(const ScatterPlotData& rhs){
		if( (this->z > rhs.z) )
		{ return 0;}
		else {return 1;}
	}
	
	inline bool operator<(const ScatterPlotData& rhs){
		if( (this->z < rhs.z) )
		{ return 0;}
		else {return 1;}
	}
};

class ScatteringModel
{
	
public:

	// Constructor
	ScatteringModel();

// Collimation Functions	
	// Set ScatterType
	void SetScatterType(int st);
	
	// Calculate the particle path length in given material using scattering processes
	double PathLength(Material* mat, double E0);
	
	// Energy loss via ionisation
	void EnergyLoss(PSvector& p, double x, Material* mat, double E0, double E1);
	
	// Advanced energy loss via ionisation
	void EnergyLoss(PSvector& p, double x, Material* mat, double E0);

	// Multiple Coulomb scattering
	void Straggle(PSvector& p, double x, Material* mat, double E1, double E2);
	
	// Function performs scattering and returns true if inelastic scatter
	bool ParticleScatter(PSvector& p, Material* mat, double E);

	// Used for output
	void DeathReport(PSvector& p, double x, double position, vector<double>& lost);

// Other Functions

	// Add/clear ScatteringProcesses
	void AddProcess(Collimation::ScatteringProcess* S){ Processes.push_back(S); fraction.push_back(0); }
	void ClearProcesses(){Processes.clear();}

	// Scatter plot
	void ScatterPlot(ParticleTracking::Particle& p, double z, int turn, string name);
	void SetScatterPlot(string name, int single_turn = 0);
	void OutputScatterPlot(string directory, int seed = 0);
	vector<string> ScatterPlotNames;
	bool ScatterPlot_on;
	vector <ScatterPlotData*> StoredScatterPlotData;
	
	// Jaw impact
	void JawImpact(ParticleTracking::Particle& p, int turn, string name);
	void SetJawImpact(string name, int single_turn = 0);
	void OutputJawImpact(string directory, int seed = 0);
	vector<string> JawImpactNames;
	bool JawImpact_on;
	vector <JawImpactData*> StoredJawImpactData;

	// vector holding all scattering processes
	vector <Collimation::ScatteringProcess*> Processes;
	// vector with fractions of the total scattering cross section assigned to each ScatteringProcess
	vector <double> fraction;
	
	//Store calculated CrossSections data to save time
	std::map< string, Collimation::CrossSections* > stored_cross_sections;
	std::map< string, Collimation::CrossSections* >::iterator CS_iterator;	
	
	int GetScatteringPhysicsModel(){return ScatteringPhysicsModel;}
		
protected:

private:

	//0 = SixTrack, 1 = ST+Ad Ion, 2 = ST + Ad El, 3 = ST + Ad SD, 4 = MERLIN	
    int ScatteringPhysicsModel;
};

} //end namespace Collimation
#endif
