/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ScatteringModel_h
#define ScatteringModel_h 1

#include <iostream>
#include <cmath>
#include <map>

#include "merlin_config.h"

#include "PSvector.h"

#include "ParticleBunch.h"

#include "Material.h"
#include "ScatteringProcess.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

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
class ScatteringModel
{

public:

	/**
	 * Constructor
	 */
	ScatteringModel();
	virtual ~ScatteringModel();

	/**
	 * Collimation Functions
	 * Set ScatterType
	 */
	void SetScatterType(int st);

	/**
	 * Calculate the particle path length in given material using scattering processes
	 */
	double PathLength(Material* mat, double E0);

	/**
	 * Dispatches to EnergyLossSimple or EnergyLossFull
	 */
	void EnergyLoss(PSvector& p, double x, Material* mat, double E0);

	/**
	 * Multiple Coulomb scattering
	 */
	void Straggle(PSvector& p, double x, Material* mat, double E1, double E2);

	/**
	 * Function performs scattering and returns true if inelastic scatter
	 */
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

protected:
	/**
	 * vector holding all scattering processes
	 */
	std::vector<Collimation::ScatteringProcess*> Processes;

	/**
	 * vector with fractions of the total scattering cross section assigned to each ScatteringProcess
	 */
	std::vector<double> fraction;

	/**
	 * Store calculated CrossSections data to save time
	 */
	std::map<std::string, Collimation::CrossSections*> stored_cross_sections;
	std::map<std::string, Collimation::CrossSections*>::iterator CS_iterator;
	EnergyLossMode energy_loss_mode;

private:

	/**
	 * Energy loss via ionisation
	 */
	void EnergyLossSimple(PSvector& p, double x, Material* mat, double E0);

	/**
	 * Advanced energy loss via ionisation
	 */

	void EnergyLossFull(PSvector& p, double x, Material* mat, double E0);
	//0 = SixTrack, 1 = ST+Ad Ion, 2 = ST + Ad El, 3 = ST + Ad SD, 4 = MERLIN
	int ScatteringPhysicsModel; // Still required for CrossSections
};

} //end namespace Collimation
#endif
