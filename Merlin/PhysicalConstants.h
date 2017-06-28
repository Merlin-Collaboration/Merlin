// Specification for PhysicalConstants
// last modified: 01/30/98 at 15:26:26
// This file is part of the MERLIN class library
// version 1.0beta
// ---------------------------------------------

// (c) 1997 Nicholas J. Walker (DESY) -- All Rights Reserved --
// ------------------------------------------------------------


#ifndef PhysicalConstants_h
#define PhysicalConstants_h 1

#include "merlin_config.h"

#include "PhysicalUnits.h"
#include "NumericalConstants.h"

namespace PhysicalConstants
{
// Data Members for Class Attributes

extern const double Avogadro;

extern const double AtomicMassUnit;

extern const double SpeedOfLight;

extern const double ElectronMass;

extern const double ProtonMass;

extern const double PlanckConstant;

extern const double PlanckConstantBar;

extern const double ElectronCharge;

extern const double ElectronMassMeV;

extern const double ProtonMassMeV;

extern const double ProtonMassGeV;

extern const double FreeSpacePermeability;

extern const double FreeSpacePermittivity;

extern const double ElectronRadius;

extern const double ElectronGe;

extern const double MuonMass;

extern const double MuonMassMeV;

extern const double MuonLifetime;

extern const double FineStructureConstant;

extern const double PionZeroMassMeV;

extern double LorentzBeta (double gamma);

extern double LorentzGamma (double beta);
extern double LorentzGamma (double momentum, double mass);
}

// Class Utility PhysicalConstants
#endif
