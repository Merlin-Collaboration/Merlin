/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: (c) 1997 Nicholas J. Walker (DESY) -- All Rights Reserved --
 */

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

extern double LorentzBeta(double gamma);

extern double LorentzGamma(double beta);
extern double LorentzGamma(double momentum, double mass);
}

// Class Utility PhysicalConstants
#endif
