// Body for PhysicalConstants
// last modified: 01/30/98 at 15:26:27
// This file is part of the MERLIN class library
// version 1.0beta
// ---------------------------------------------

// (c) 1997 Nicholas J. Walker (DESY) -- All Rights Reserved --
// ------------------------------------------------------------


#include <cmath>

#include "NumericalUtils/PhysicalUnits.h"

// PhysicalConstants
#include "NumericalUtils/PhysicalConstants.h"


// Class Utility PhysicalConstants
namespace PhysicalConstants {
using namespace PhysicalUnits;


const double SpeedOfLight = 2.997925e+08*meter/second;



const double ElectronMass = 9.10956e-31; // kilogram



const double ProtonMass = 1.67261e-27; // kilogram



const double PlanckConstant = 6.62620e-34; // Joule second



const double ElectronCharge = 1.60219e-19; // Coulomb



const double ElectronMassMeV = 0.511005; // MeV



const double ProtonMassMeV = 938.2578; // MeV


const double FreeSpacePermeability = 16.0e-7*atan(1.0); // Henry per meter

const double FreeSpacePermittivity = 1.0/FreeSpacePermeability/SpeedOfLight/SpeedOfLight; // Farad per meter

const double ElectronRadius = ElectronCharge*ElectronCharge/16.0/atan(1.0)/FreeSpacePermittivity/ElectronMass/SpeedOfLight/SpeedOfLight; // meter

const double ElectronGe = 0.00115965219;

};

