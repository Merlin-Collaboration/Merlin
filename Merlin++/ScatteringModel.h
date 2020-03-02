/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ScatteringModel_h
#define ScatteringModel_h 1

#include <iostream>
#include <string>
#include <map>

#include "merlin_config.h"
#include "PSvector.h"
#include "ScatteringProcess.h"
#include "utils.h"

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
		if((this->z != rhs.z))
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
		if((this->z > rhs.z))
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
		if((this->z < rhs.z))
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

};

enum EnergyLossMode
{
	SimpleEnergyLoss,
	FullEnergyLoss

};


/**
 * Base class for scattering models
 *
 * The user can customise a ScatteringModel using AddProcess(), or can
 * use the predefined models such as ScatteringModelMerlin.
 */

struct ScatterModelDetails
{
	std::vector<double> Xsection = std::vector<double>(5, 0.0);
	std::vector<Collimation::ScatteringProcess*> Processes{0, 0, 0, 0, 0, 0};

};

class ScatteringModel
{

public:
	MaterialProperties* oldMaterial;   // keep track of changing collimators
        int ModelType; // 0 is 'Merlin', 1 is 'Sixtrack'.  Not elegant !! FIX!!
	map<MaterialProperties*, ScatterModelDetails*> saveDetails;
	/**
	 * Constructor
	 */
	ScatteringModel(int model=0);
	virtual ~ScatteringModel();
	void Configure(MaterialProperties *, double Energy);    // material not known at
	// construct time and may change

	/**
	 * Collimation Functions
	 * Set ScatterType
	 */
	void SetScatterType(int st);

	/**
	 * Calculate the particle path length in given material using scattering processes
	 */
	double PathLength(MaterialProperties* mat, double E0);

	/**
	 * Dispatches to EnergyLossSimple or EnergyLossFull
	 */
	void EnergyLoss(PSvector& p, double x, MaterialProperties* mat, double E0);

	/**
	 * Multiple Coulomb scattering
	 */
	void Straggle(PSvector& p, double x, MaterialProperties* mat, double E1, double E2);

	/**
	 * Function performs scattering and returns true if inelastic scatter
	 */
	bool ParticleScatter(PSvector& p, MaterialProperties* mat, double E);

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
	std::vector<ScatterPlotData*> StoredScatterPlotData;

	// Jaw impact
	void JawImpact(ParticleTracking::Particle& p, int turn, std::string name);
	void SetJawImpact(std::string name, int single_turn = 0);
	void OutputJawImpact(std::string directory, int seed = 0);
	std::vector<std::string> JawImpactNames;
	bool JawImpact_on;
	std::vector<JawImpactData*> StoredJawImpactData;

	int GetScatteringPhysicsModel()
	{
		return ScatteringPhysicsModel;
	}

public:
	/**
	 * vector holding all scattering processes
	 */
	// adapted RJB	std::vector<Collimation::ScatteringProcess*> Processes;
	std::vector<Collimation::ScatteringProcess*> Processes{0, 0, 0, 0, 0 ,0};

	/**
	 * vector with fractions of the total scattering cross section assigned to each ScatteringProcess
	 */
	std::vector<double> fraction;
// RJB   To supersede fraction
	std::vector<double> Xsection = std::vector<double>(5, 0.0);

	EnergyLossMode energy_loss_mode;

private:

	/**
	 * Energy loss via ionisation
	 */
	void EnergyLossSimple(PSvector& p, double x, MaterialProperties* mat, double E0);

	/**
	 * Advanced energy loss via ionisation
	 */

	void EnergyLossFull(PSvector& p, double x, MaterialProperties* mat, double E0);
	//0 = SixTrack, 1 = ST+Ad Ion, 2 = ST + Ad El, 3 = ST + Ad SD, 4 = MERLIN
	int ScatteringPhysicsModel; // Still required for CrossSections
};

} //end namespace Collimation
#endif
