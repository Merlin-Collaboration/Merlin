#ifndef ScatteringModel_h
#define ScatteringModel_h 1

#include <iostream>
#include <cmath>
#include <map>

#include "merlin_config.h"

#include "BeamModel/PSvector.h"

#include "Collimators/Material.h"
#include "Collimators/ScatteringProcess.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"


namespace Collimation {

class ScatteringModel
{
	
public:

	// Constructor
	ScatteringModel();

	// Add ScatteringProcesses
	void AddProcess(Collimation::ScatteringProcess* S){ Processes.push_back(S); fraction.push_back(0); }
	void ClearProcesses(){Processes.clear();}
	
	// Set ScatterType
	void SetScatterType(int st);
	
	// Calculate the particle path length in given material using scattering processes
	double PathLength(Material* mat, double E0);

	// Function performs scattering and returns true if inelastic scatter
	bool ParticleScatter(PSvector& p, Material* mat, double E0);
	//~ bool ParticleScatter(PSvector& p, Material* mat, double E0, double sigpNtot, double sigR);

	// Used for output
	void DeathReport(PSvector& p, double x, double position, vector<double>& lost);

	//void Straggle(PSvector& p, double x, MaterialProperties* prop, double E0);
	//void EnergyLoss(PSvector& p, double x, MaterialProperties* prop, double E0);
	
	// Energy loss via ionisation
	void EnergyLoss(PSvector& p, double x, Material* mat, double E0, double E1);
	
	// Advanced energy loss via ionisation
	void EnergyLoss(PSvector& p, double x, Material* mat, double E0);

	// Multiple Coulomb scattering
	void Straggle(PSvector& p, double x, Material* mat, double E1, double E2);

	// vector holding all scattering processes
	vector <Collimation::ScatteringProcess*> Processes;
	// vector wil fractions of the total scattering cross section assigned to each ScatteringProcess
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
