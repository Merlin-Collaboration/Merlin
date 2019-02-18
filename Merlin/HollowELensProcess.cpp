/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>

#include "utils.h"
#include "RandomNG.h"
#include "HollowELensProcess.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace std;

namespace ParticleTracking
{

HollowELensProcess::HollowELensProcess(int priority) :
	ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), currentComponentHEL(nullptr)
{
}

void HollowELensProcess::InitialiseProcess(Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	if(!currentBunch)
	{
		cout << "HEL warning: !currentBunch" << endl;
	}
}

void HollowELensProcess::SetCurrentComponent(AcceleratorComponent& component)
{
	HollowElectronLens* aHollowELens = dynamic_cast<HollowElectronLens*>(&component);
	active = (currentBunch != nullptr) && (aHollowELens);

	if(active)
	{
		currentComponent = &component;
		currentComponentHEL = aHollowELens;

		double Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
	}
	else
	{
		currentComponent = nullptr;
	}
}

void HollowELensProcess::DoProcess(double /*ds*/)
{
	// NB
	// CalcThetaMax returns +ve theta
	// CalcKick Radial and Simple return -ve theta
	// DoProcess should have p.xp() += theta * sin(ParticleAngle)

	double theta = 0;
	double Gamma_p = 0;
	double ParticleAngle;

	bool SimpleProfile = currentComponentHEL->SimpleProfile;
	bool ACSet = currentComponentHEL->ACSet;

	ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(), currentBunch->GetTotalCharge()
		/ currentBunch->size());
	newbunch->clear();
	newbunch->swap(*currentBunch);

	if(ProtonBeta == 0)
	{
		Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
	}

	// Have to increment Turn as the process doesn't have access to the turn value from user code
	currentComponentHEL->Turn++;

	switch(currentComponentHEL->OMode)
	{
	case DC:
	{
		//HEL always on
		for(PSvectorArray::iterator p = newbunch->begin(); p != newbunch->end(); p++)
		{
			if(SimpleProfile)
			{
				theta = CalcKickSimple(*p);
			}
			else
			{
				theta = CalcKickRadial(*p);
			}

			if(theta != 0)
			{
				ParticleAngle = atan2((*p).y(), (*p).x());
				//~ // Particle phase space angle and amplitude (radius)
				(*p).xp() += theta * cos(ParticleAngle);
				(*p).yp() += theta * sin(ParticleAngle);
			}
		}
	}
	break;
	case AC:
	{
		// Resonant HEL kick - Adapted from V. Previtali's SixTrack elense
		if(ACSet)
		{
			double TuneVarPerStep = currentComponentHEL->TuneVarPerStep;
			double DeltaTune = currentComponentHEL->DeltaTune;
			double MinTune = currentComponentHEL->MinTune;
			int Turn = currentComponentHEL->Turn;
			double TurnsPerStep = currentComponentHEL->TurnsPerStep;
			double Nstep = currentComponentHEL->Nstep;
			double Tune = currentComponentHEL->Tune;
			double Multiplier = currentComponentHEL->Multiplier;
			for(PSvectorArray::iterator p = newbunch->begin(); p != newbunch->end(); p++)
			{
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}

				if(theta != 0)
				{
					double OpTune;
					if((TuneVarPerStep != 0) && (DeltaTune != 0))
					{
						OpTune = MinTune + fmod((floor(Turn / TurnsPerStep)), (Nstep)) * TuneVarPerStep;
					}
					else
					{
						OpTune = Tune;
					}

					double Phi = Multiplier * (Turn * OpTune * 2 * pi);
					theta *= 0.5 * (1 + cos(Phi));

					ParticleAngle = atan2((*p).y(), (*p).x());
					// Particle phase space angle and amplitude (radius)
					(*p).xp() += theta * cos(ParticleAngle);
					(*p).yp() += theta * sin(ParticleAngle);
				}
			}
		}
		else
		{
			cout << "\n\tHEL Warning: AC variables not set" << endl;
		}
	}
	break;
	case Diffusive:
	{
		// HEL randomly switched on/off on a turn by turn basis
		double rando = RandomNG::uniform(-1, 1);

		if(rando >= 0)
		{
			for(PSvectorArray::iterator p = newbunch->begin(); p != newbunch->end(); p++)
			{
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}

				ParticleAngle = atan2((*p).y(), (*p).x());
				if(theta != 0)
				{
					// Particle phase space angle and amplitude (radius)
					(*p).xp() += theta * cos(ParticleAngle);
					(*p).yp() += theta * sin(ParticleAngle);
				}
			}
		}
	}
	break;
	case Turnskip:
	{
		int SkipTurn = currentComponentHEL->SkipTurn;
		int Turn = currentComponentHEL->Turn;
		// HEL switched on/off if turn = muliple of n
		if(SkipTurn == 0)
		{
			cout << "\n\tHEL warning: SkipTurn not set, autoset to 2" << endl;
			SkipTurn = 2;
		}
		if((Turn % SkipTurn) == 0)
		{
			for(PSvectorArray::iterator p = newbunch->begin(); p != newbunch->end(); p++)
			{
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}

				ParticleAngle = atan2((*p).y(), (*p).x());
				if(theta != 0)
				{
					// Particle phase space angle and amplitude (radius)
					(*p).xp() += theta * cos(ParticleAngle);
					(*p).yp() += theta * sin(ParticleAngle);
				}
			}
		}
	}
	break;
	} //end switch

	//Apply changes to particle bunch
	currentBunch->clear();
	currentBunch->swap(*newbunch);
	delete newbunch;
}

double HollowELensProcess::GetMaxAllowedStepSize() const
{
	return currentComponent->GetLength();
}

double HollowELensProcess::CalcThetaMax(double r)
{
	if(r == 0)
	{
		return 0;
	}

	bool ElectronDirection = currentComponentHEL->ElectronDirection;
	double EffectiveLength = currentComponentHEL->EffectiveLength;
	double Current = currentComponentHEL->Current;
	double ElectronBeta = currentComponentHEL->ElectronBeta;
	double Rigidity = currentComponentHEL->Rigidity;

	double ThetaMax;
	if(ElectronDirection)
	{
		// HEL electrons travelling opposite to LHC protons (summed kick)
		ThetaMax = (2 * EffectiveLength * Current * (1 + (ElectronBeta * ProtonBeta))) / (r * 1E7 * Rigidity
			* ElectronBeta * ProtonBeta);
	}
	else
	{
		// HEL electrons travelling in the same direction to LHC protons (smaller kick and opposite)
		ThetaMax = -(2 * EffectiveLength * Current * (1 - (ElectronBeta * ProtonBeta))) / (r * 1E7 * Rigidity
			* ElectronBeta * ProtonBeta);
	}

	return ThetaMax;
}

double HollowELensProcess::CalcKickSimple(Particle &p)
{
	// Start of HEL
	double x = p.x();
	double y = p.y();

	double XOffset = currentComponentHEL->XOffset;
	double YOffset = currentComponentHEL->YOffset;

	// Calculate particle transverse vector ('radius' in xy space)
	double R = sqrt(pow((x - XOffset), 2) + pow((y - YOffset), 2));
	return CalcKickSimple(R);
}

double HollowELensProcess::CalcKickSimple(double R)
{
	double thet = 0;
	double Rmin = currentComponentHEL->GetRmin();
	double Rmax = currentComponentHEL->GetRmax();

	if(R <= Rmin)
	{
		return 0;
	}
	else if(R < Rmax && R > Rmin)
	{
		thet = ((pow(R, 2) - pow(Rmin, 2)) / (pow(Rmax, 2) - pow(Rmin, 2))) * CalcThetaMax(R);
		return -thet;
	}
	else
	{
		thet = CalcThetaMax(R);
		return -thet;
	}
}

double HollowELensProcess::CalcKickRadial(Particle &p)
{
	// Start of HEL
	double x = p.x();
	double y = p.y();

	double XOffset = currentComponentHEL->XOffset;
	double YOffset = currentComponentHEL->YOffset;

	// Calculate particle transverse vector ('radius' in xy space)
	double R = sqrt(pow((x - XOffset), 2) + pow((y - YOffset), 2));
	return CalcKickRadial(R);
}

double HollowELensProcess::CalcKickRadial(double R)
{
	double f = 0;
	double thet = 0;
	double Rmin = currentComponentHEL->GetRmin();
	bool LHC_Radial = currentComponentHEL->LHC_Radial;

	// Adapted from V. Previtali's SixTrack elense implementation
	if(R <= Rmin)
	{
		return 0;
	}

	// Define boundaries between parameterisation of measured radial profile
	double r0, r1, r2, r3, r4;

	if(!LHC_Radial)
	{
		// Tevatron HEL 1.2A, 2m, 5KeV, 4-6.8sig
		r0 = 222.5;
		r1 = 252.5;
		r2 = 287;
		r3 = 364.5;
		r4 = 426.5;
	}
	else
	{
		// LHC HEL 5A, 3m, 10KeV, 4-8sig
		r1 = 265;       // 0 - initial rise				(x1.191)
		r2 = 315;       // rise - straight section		(x1.248)
		r3 = 435;       // straight - left of peak		(x1.193)
		r4 = 505;       // left - right of peak			(x1.184)
	}

	double elense_r_min = Rmin; //Need to calculate 4 sigma at this point

	double x0 = elense_r_min;
	double y0 = 0;
	double x1 = r1 / r0 * elense_r_min;
	double y1 = 917;
	double x2 = r2 / r0 * elense_r_min;
	double y2 = 397;
	double x3 = r3 / r0 * elense_r_min;
	double y3 = 228;
	double x4 = r4 / r0 * elense_r_min;
	double y4 = 0;

	double n0 = ((y1 - y0) / (x1 - x0) / 3) * pow(x1, 3) + (y0 - x0 * (y1 - y0) / (x1 - x0)) * pow(x1, 2) / 2 - (((y1
		- y0) / (x1 - x0) / 3) * pow(x0, 3) + (y0 - x0 * (y1 - y0) / (x1 - x0)) * pow(x0, 2) / 2);
	double n1 = (y2 - y1) / (x2 - x1) * pow(x2, 3) / 3 + (y1 - x1 * (y2 - y1) / (x2 - x1)) * pow(x2, 2) / 2 - ((y2
		- y1) / (x2 - x1) * pow(x1, 3) / 3 + (y1 - x1 * (y2 - y1) / (x2 - x1)) * pow(x1, 2) / 2);
	double n2 = (y3 - y2) / (x3 - x2) * pow(x3, 3) / 3 + (y2 - x2 * (y3 - y2) / (x3 - x2)) * pow(x3, 2) / 2 - ((y3
		- y2) / (x3 - x2) * pow(x2, 3) / 3 + (y2 - x2 * (y3 - y2) / (x3 - x2)) * pow(x2, 2) / 2);
	double n3 = (y4 - y3) / (x4 - x3) * pow(x4, 3) / 3 + (y3 - x3 * (y4 - y3) / (x4 - x3)) * pow(x4, 2) / 2 - ((y4
		- y3) / (x4 - x3) * pow(x3, 3) / 3 + (y3 - x3 * (y4 - y3) / (x4 - x3)) * pow(x3, 2) / 2);
	double ntot = n0 + n1 + n2 + n3;

	if(R < x0)
	{
		f = 0;
		cout << "HEL warning: Radial profile: R < x0" << endl;
	}
	else if(R < x1)
	{
		f = (((y1 - y0) / (x1 - x0) / 3) * pow(R, 3) + (y0 - x0 * (y1 - y0) / (x1 - x0)) * pow(R, 2) / 2 - (((y1 - y0)
			/ (x1 - x0) / 3) * pow(x0, 3) + (y0 - x0 * (y1 - y0) / (x1 - x0)) * pow(x0, 2) / 2)) / ntot;
	}
	else if(R < x2)
	{
		f = (n0 + (y2 - y1) / (x2 - x1) * pow(R, 3) / 3 + (y1 - x1 * (y2 - y1) / (x2 - x1)) * pow(R, 2) / 2 - ((y2
			- y1) / (x2 - x1) * pow(x1, 3) / 3 + (y1 - x1 * (y2 - y1) / (x2 - x1)) * pow(x1, 2) / 2)) / ntot;
	}
	else if(R < x3)
	{
		f = (n0 + n1 + (y3 - y2) / (x3 - x2) * pow(R, 3) / 3 + (y2 - x2 * (y3 - y2) / (x3 - x2)) * pow(R, 2) / 2
			- ((y3 - y2) / (x3 - x2) * pow(x2, 3) / 3 + (y2 - x2 * (y3 - y2) / (x3 - x2)) * pow(x2, 2) / 2)) / ntot;
	}
	else if(R < x4)
	{
		f = (n0 + n1 + n2 + (y4 - y3) / (x4 - x3) * pow(R, 3) / 3 + (y3 - x3 * (y4 - y3) / (x4 - x3)) * pow(R, 2) / 2
			- ((y4 - y3) / (x4 - x3) * pow(x3, 3) / 3 + (y3 - x3 * (y4 - y3) / (x4 - x3)) * pow(x3, 2) / 2)) / ntot;
	}
	else
	{
		f = 1;
	}

	thet = -CalcThetaMax(R) * f;
	return thet;
}

// TODO: Should reimplement a HollowElectionLensConfiguration system (like ApertureConfiguration)

void HollowELensProcess::OutputProfile(std::ostream* os, double E, double min, double max)
{
	// Have to set ProtonBeta if running before it is set
	double Gamma_p = LorentzGamma(E, ProtonMass);
	ProtonBeta = LorentzBeta(Gamma_p);

	if(!currentComponentHEL)
	{
		cerr << "currentComponentHEL not set. Call SetCurrentComponent()" << endl;
		return;
	}
	double Rmin = currentComponentHEL->GetRmin();
	double Rmax = currentComponentHEL->GetRmax();
	double Sigma_x = currentComponentHEL->Sigma_x;
	double EffectiveLength = currentComponentHEL->EffectiveLength;
	double Current = currentComponentHEL->Current;
	double Rigidity = currentComponentHEL->Rigidity;
	double ElectronBeta = currentComponentHEL->ElectronBeta;

	cout << " Rmin = " << Rmin << ", = " << Rmin / Sigma_x << " sigma" << endl;
	cout << " Rmax = " << Rmax << ", = " << Rmax / Sigma_x << " sigma" << endl;
	cout << " L = " << EffectiveLength << endl;
	cout << " Current = " << Current << endl;
	cout << " Brho = " << Rigidity << endl;
	cout << " ElectronBeta = " << ElectronBeta << endl;
	cout << " ProtonBeta = " << ProtonBeta << endl;

	double r = 0;
	int points = 1000;

	(*os) << "#r\tkick_radial\tkick_simple\t|kick_r|\t|kick_s|" << endl;

	for(int i = 1; i < points; ++i)
	{
		(*os) << (r / Sigma_x) << "\t" << CalcKickRadial(r) << "\t" << CalcKickSimple(r) << "\t" << sqrt(pow(
				CalcKickRadial(r), 2)) << "\t" << sqrt(pow(CalcKickSimple(r), 2)) << endl;
		r += ((max - min) / points) * Sigma_x;
	}
}

}
