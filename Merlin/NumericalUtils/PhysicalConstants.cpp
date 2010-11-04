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

//Updated with PDG 2008 values

const double Avogadro = 6.02214179e23;

const double AtomicMassUnit = 931.494028*MeV; //MeV given by PDG 2008

const double SpeedOfLight = 2.99792458e+08*meter/second;

const double ElectronMass = 9.10938215e-31; // kilogram (old = 9.10956e-31)

const double ProtonMass = 1.672621637e-27; // kilogram  (old = 1.67261e-27)

const double PlanckConstant = 6.62606896e-34; // Joule second (old = 6.62620e-34)

const double ElectronCharge = 1.602176487e-19; // Coulomb (old = 1.60219e-19)

const double ElectronMassMeV = 0.510998910; // MeV (old = 0.511005)

const double ProtonMassMeV = 938.272013; // MeV (old = 938.2578)

const double FreeSpacePermeability = 16.0e-7*atan(1.0); // Henry per meter

const double FreeSpacePermittivity = 1.0/FreeSpacePermeability/SpeedOfLight/SpeedOfLight; // Farad per meter

const double ElectronRadius = ElectronCharge*ElectronCharge/16.0/atan(1.0)/FreeSpacePermittivity/ElectronMass/SpeedOfLight/SpeedOfLight; // meter

const double ElectronGe = 0.001159652181; //Magnetic moment anomaly e cm (old = 0.00115965219)

const double MuonMass = ProtonMass * MuonMassMeV/ProtonMassMeV; //AtomicMassUnit * 0.11342892;

const double MuonMassMeV = 105.65836;

const double MuonLifetime = 2.1970e-6 * second;	//In the rest frame

const double FineStructureConstant = (ElectronCharge * ElectronCharge * SpeedOfLight * FreeSpacePermeability) / (2 * PlanckConstant);
}
