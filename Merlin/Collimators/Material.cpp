#include "Collimators/Material.hpp"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

double material::CalculateElectronDensity()
{
	//rho is g/cm^3
	return atomic_number * Avogadro * rho / (A * pow(centimeter,3)); // n_e m^-3 (1e6 conversion from cm^3)
}

double material::CalculatePlasmaEnergy()
{
	//returns the Plasma energy in GeV
	return (PlanckConstantBar * sqrt((ElectronDensity * pow(ElectronCharge,2)) / (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
}


double material::CalculateMeanExcitationEnergy()
{
	return atomic_number * 10.0 *eV;
}
