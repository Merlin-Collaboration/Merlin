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

namespace Collimation
{

struct JawImpactData
{
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

struct ScatterPlotData
{
	int turn;
	int ID;
	double x;
	double xp;
	double y;
	double yp;
	double z;
	string name;

	inline bool operator==(const ScatterPlotData& rhs)
	{
		if( (this->z != rhs.z) )
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

	inline bool operator>(const ScatterPlotData& rhs)
	{
		if( (this->z > rhs.z) )
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

	inline bool operator<(const ScatterPlotData& rhs)
	{
		if( (this->z < rhs.z) )
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}
};

/**
 * Base class for scattering models
 *
 * The user can customise a ScatteringModel using AddProcess(), or can
 * use the predefiend models such as ScatteringModelMerlin.
 */
class ScatteringModel
{

public:

	// Constructor
	ScatteringModel();
	virtual ~ScatteringModel();

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

// Other Functions

	/**
	 * Add ScatteringProcesses to the model
	 */
	virtual void AddProcess(Collimation::ScatteringProcess* S)
	{
		Processes.push_back(S);
		fraction.push_back(0);
	}
	void ClearProcesses()
	{
		Processes.clear();
	}

	// Scatter plot
	void ScatterPlot(ParticleTracking::Particle& p, double z, int turn, std::string name);
	void SetScatterPlot(std::string name, int single_turn = 0);
	void OutputScatterPlot(std::string directory, int seed = 0);
	std::vector<std::string> ScatterPlotNames;
	bool ScatterPlot_on;
	std::vector <ScatterPlotData*> StoredScatterPlotData;

	// Jaw impact
	void JawImpact(ParticleTracking::Particle& p, int turn, std::string name);
	void SetJawImpact(std::string name, int single_turn = 0);
	void OutputJawImpact(std::string directory, int seed = 0);
	std::vector<std::string> JawImpactNames;
	bool JawImpact_on;
	std::vector <JawImpactData*> StoredJawImpactData;

	int GetScatteringPhysicsModel()
	{
		return ScatteringPhysicsModel;
	}

protected:
	// vector holding all scattering processes
	std::vector <Collimation::ScatteringProcess*> Processes;

	// vector with fractions of the total scattering cross section assigned to each ScatteringProcess
	std::vector <double> fraction;

	//Store calculated CrossSections data to save time
	std::map< std::string, Collimation::CrossSections* > stored_cross_sections;
	std::map< std::string, Collimation::CrossSections* >::iterator CS_iterator;

private:

	//0 = SixTrack, 1 = ST+Ad Ion, 2 = ST + Ad El, 3 = ST + Ad SD, 4 = MERLIN
	int ScatteringPhysicsModel;
};

} //end namespace Collimation
#endif
