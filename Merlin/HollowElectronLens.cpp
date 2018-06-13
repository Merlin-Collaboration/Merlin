/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "HollowElectronLens.h"
#include "ComponentTracker.h"

const int HollowElectronLens::ID = UniqueIndex();

HollowElectronLens::HollowElectronLens(const string& id, double len, int mode, double current, double beta_e, double
	rigidity, double length_e) :
	Drift(id, len), Current(current), ElectronBeta(beta_e), Rigidity(rigidity), EffectiveLength(length_e), XOffset(0),
	YOffset(0), Turn(0), SkipTurn(0), ACSet(0), SimpleProfile(1), AlignedToOrbit(0), ElectronDirection(1), LHC_Radial(
		0)
{
	if(mode == 0)
	{
		OMode = DC;
	}
	else if(mode == 1)
	{
		OMode = AC;
	}
	else if(mode == 2)
	{
		OMode = Diffusive;
	}
	else if(mode == 3)
	{
		OMode = Turnskip;
	}
	else
	{
		cout
			<<
			"\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip"
			<< endl;
	}
}

const string& HollowElectronLens::GetType() const
{
	_TYPESTR(HollowElectronLens);
}

int HollowElectronLens::GetIndex() const
{
	return ID;
}

void HollowElectronLens::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, Drift); //HR change to stop collimator tracking 29.10.13
	//_PREPTRACK(aTracker,AcceleratorComponent);
}

void HollowElectronLens::RotateY180()
{
	// nothing to do
}

void HollowElectronLens::SetRmax(double rmax)
{
	Rmax = rmax;
}

void HollowElectronLens::SetRmin(double rmin)
{
	Rmin = rmin;
}

void HollowElectronLens::SetRadii(double rmin, double rmax)
{
	SetRmin(rmin);
	SetRmax(rmax);
}

ModelElement* HollowElectronLens::Copy() const
{
	return new HollowElectronLens(*this);
}

void HollowElectronLens::SetAC(double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi)
{
	Tune = tune;
	DeltaTune = deltatune;
	TuneVarPerStep = tunevarperstep;
	TurnsPerStep = turnsperstep;
	Multiplier = multi;
	MinTune = Tune - DeltaTune;
	Nstep = (2 * DeltaTune / TuneVarPerStep) + 1;
	Turn = 0;
	ACSet = 1;
	OMode = AC;
}

void HollowElectronLens::SetTurnskip(int skip)
{
	SkipTurn = skip;
	OMode = Turnskip;
}

void HollowElectronLens::SetElectronDirection(bool dir)
{
	ElectronDirection = dir;
	if(ElectronDirection)
	{
		cout << "HELProcess: electrons travelling opposite to protons: negative (focussing) kick" << endl;
	}
	else
	{
		cout << "HELProcess: electrons travelling in the same direction as protons: positive (defocussing) kick"
			 << endl;
	}
}
