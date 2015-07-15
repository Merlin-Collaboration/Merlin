#ifndef ScatteringModel_h
#define ScatteringModel_h 1

#include <iostream>

#include "merlin_config.h"

#include "BeamModel/PSvector.h"

#include "Collimators/Material.h"
#include "Collimators/ScatteringProcess.h"

class ScatteringModel
{
	
public:

	// Constructor
	ScatteringModel();

	// Add ScatteringProcesses
	void AddProcess(ScatteringProcess* S){ Processes.push_back(S); fraction.push_back(0); }
	
	// Calculate the particle path length in given material using scattering processes
	double PathLength(Material* mat, double E0);

	// Function performs scattering and returns true if inelastic scatter
	bool ParticleScatter(PSvector& p, Material* mat, double E0);

	// Used for output
	void DeathReport(PSvector& p, double x, double position, vector<double>& lost);

	//void Straggle(PSvector& p, double x, MaterialProperties* prop, double E0);
	//void EnergyLoss(PSvector& p, double x, MaterialProperties* prop, double E0);
	
	// Energy loss via ionisation
	double EnergyLoss(PSvector& p, double x, Material* mat, double E0, double E1);

	// Multiple Coulomb scattering
	void Straggle(PSvector& p, double x, Material* mat, double E0, double E2);

	// vector holding all scattering processes
	vector <ScatteringProcess*> Processes;
	// vector wil fractions of the total scattering cross section assigned to each ScatteringProcess
	vector <double> fraction;
	// MaterialProperties of the previous collimator (no need to reinitialise all values)
	Material* lastmat;


protected:

private:


};
#endif
