/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 5.01 (2015)
//
// Copyright: see Merlin/copyright.txt
//
// Created:		2014 HR
// Modified:
// Last Edited: 15.12.15 HR
//
/////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "utils.h"

#include "RandomNG.h"

#include "HollowELensProcess.h"

#include "LatticeFunctions.h"
#include "PhaseAdvance.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace std;

namespace ParticleTracking
{


HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity),\
	  ACSet(0), SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0), ElectronDirection(1),\
	  LHC_Radial(0), YShift(0.), XShift(0.), Elliptical(0), EllipticalSet(0), g(0.)
{
	if (mode == 0)
	{
		OMode = DC;
	}
	else if (mode == 1)
	{
		OMode = AC;
	}
	else if (mode == 2)
	{
		OMode = Diffusive;
	}
	else if (mode == 3)
	{
		OMode = Turnskip;
	}
	else
	{
		cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;
	}
}

HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity),\
	  EffectiveLength(length_e), ACSet(0), SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0),\
	  ElectronDirection(1),  LHC_Radial(0), YShift(0.), XShift(0.), Elliptical(0), EllipticalSet(0), g(0.)
{
	if (mode == 0)
	{
		OMode = DC;
	}
	else if (mode == 1)
	{
		OMode = AC;
	}
	else if (mode == 2)
	{
		OMode = Diffusive;
	}
	else if (mode == 3)
	{
		OMode = Turnskip;
	}
	else
	{
		cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;
	}

}

HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double rmin,\
                                        double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity), ACSet(0),\
	  SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0), ElectronDirection(1), LHC_Radial(0),\
	  YShift(0.), XShift(0.), Elliptical(0), EllipticalSet(0), g(0.)
{
	if (mode == 0)
	{
		OMode = DC;
	}
	else if (mode == 1)
	{
		OMode = AC;
	}
	else if (mode == 2)
	{
		OMode = Diffusive;
	}
	else if (mode == 3)
	{
		OMode = Turnskip;
	}
	else
	{
		cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;
	}
	SetRadiiSigma(rmin, rmax, model, emittance_x, emittance_y, twiss);
}

void HollowELensProcess::InitialiseProcess (Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	if(!currentBunch)
	{
		cout << "HEL warning: !currentBunch" << endl;
	}
}

void HollowELensProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	HollowElectronLens* aHollowELens = dynamic_cast<HollowElectronLens*>(&component);
	active = (currentBunch!=nullptr) && (aHollowELens);

	if(active)
	{
		currentComponent = &component;
		//~ EffectiveLength = currentComponent->GetLength();

		double Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
	}
	else
	{
		currentComponent = nullptr;
	}
}

void HollowELensProcess::SetAC (double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi)
{
	Tune = tune;
	DeltaTune = deltatune;
	TuneVarPerStep = tunevarperstep;
	TurnsPerStep = turnsperstep;
	Multiplier = multi;
	MinTune = Tune - DeltaTune;
	MaxTune = Tune + DeltaTune;
	Nstep = (2 * DeltaTune / TuneVarPerStep)+1;
	Turn = 0;
	ACSet = 1;
	OMode = AC;
}

void HollowELensProcess::SetTurnskip (int skip)
{
	SkipTurn = skip;
	OMode = Turnskip;
}

void HollowELensProcess::DoProcess (double ds)
{

	// NB
	// CalcThetaMax returns +ve theta
	// CalcKick Radial and Simple return -ve theta
	// DoProcess should have p.xp() += theta * sin(ParticleAngle)

	double theta = 0;
	double Gamma_p = 0;

	ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(), currentBunch->GetTotalCharge()/currentBunch->size());
	newbunch->clear();
	newbunch->swap(*currentBunch);

	if(ProtonBeta == 0)
	{
		Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
	}

	// Have to increment Turn as the process doesn't have access to the turn value from user code
	++Turn;

	switch (OMode)
	{
	case DC:
	{
		//HEL always on
		for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++)
		{
			//~ cout << "\n\tDC Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
			if(SimpleProfile)
			{
				theta = CalcKickSimple(*p);
			}
			else
			{
				theta = CalcKickRadial(*p);
			}

			//~ cout << "\n\tX = " << (*p).x() << endl;
			//~ cout << "\n\tY = " << (*p).y() << endl;
			//~ cout << "\n\tAngle = " << ParticleAngle << endl;
			//~ cout << "\n\tRadius = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2)) << endl;
			//~ cout << "\n\tRadius = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2))/ sqrt(pow(293.031E-6,2) + pow(287.801E-6,2))<< " sigma" << endl;

			if(theta!=0)
			{

				if(!Elliptical)
				{
					ParticleAngle = atan2((*p).y(), (*p).x());
				}
				else
				{
					ParticleAngle = atan2( ((*p).y()-YShift) , ((*p).x()-XShift)  );
				}

				//~ cout << "\n Elliptical = " << Elliptical << endl;
				//~ cout << " x = " << (*p).x() << endl;
				//~ cout << " y = " << (*p).y() << endl;
				//~ cout << " xp = " << (*p).xp() << endl;
				//~ cout << " yp = " << (*p).yp() << endl;

				//~ cout << " SemiMajor = " << SemiMajor << endl;
				//~ cout << " SemiMinor = " << SemiMinor << endl;
				//~ cout << " x_shift = " << XShift << endl;
				//~ cout << " y_shift = " << YShift << endl;
				//~ cout << " x - x_shift = " << (*p).x() - XShift << endl;
				//~ cout << " y - y_shift = " << (*p).y() - YShift << endl;

				//~ cout << " ParticleAngle = " << ParticleAngle << endl;

				//~ cout << " Rmin = " << Rmin << endl;
				//~ cout << " Rmax = " << Rmax << endl;

				//~ cout << " R_beam = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2)) << endl;
				//~ cout << " R_HEL = " <<sqrt( pow(((*p).x()-XShift),2) + pow(((*p).y()-YShift),2) ) << endl;

				//~ cout << " L = " << EffectiveLength << endl;
				//~ cout << " max_kick = " << CalcThetaMax(R) << endl;
				//~ cout << " Current = " << Current << endl;
				//~ cout << " Brho = " << Rigidity << endl;
				//~ cout << " Gamma_p = " << Gamma_p << endl;
				//~ cout << " theta = " << theta << endl;
				//~ cout << "\n" << endl;

				//~ // Particle phase space angle and amplitude (radius)
				(*p).xp() += theta * cos(ParticleAngle);
				(*p).yp() += theta * sin(ParticleAngle);
			}
		}
	}
	break;
	case AC:
	{
		//~ cout << "\n\tAC Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
		// Resonant HEL kick - Adapted from V. Previtali's SixTrack elense
		if(ACSet)
		{
			for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++)
			{
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}

				if(theta!=0)
				{

					if( (TuneVarPerStep !=0) && (DeltaTune !=0) )
					{
						OpTune = MinTune + fmod((floor(Turn/TurnsPerStep)),(Nstep));
					}
					else
					{
						OpTune = Tune;
					}

					Phi = Multiplier * ( Turn * OpTune * 2 * pi );
					theta *= 0.5*(1 + cos(Phi));

					if(!Elliptical)
					{
						ParticleAngle = atan2((*p).y(), (*p).x());
					}
					else
					{
						ParticleAngle = atan2( ((*p).y()-YShift) , ((*p).x()-XShift)  );
					}

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
		//~ cout << "\n\tDiffusive Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
		// HEL randomly switched on/off on a turn by turn basis
		double rando = RandomNG::uniform(-1,1);

		if (rando >=0)
		{
			for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++)
			{
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}


				//~ if ((*p).x() < 0){	ParticleAngle = pi + atan2((*p).y(), (*p).x());	}
				//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}
				//~ if ((*p).x() < 0){	ParticleAngle = 2*pi + atan2((*p).y(), (*p).x());	}
				//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}
				if(!Elliptical)
				{
					ParticleAngle = atan2((*p).y(), (*p).x());
				}
				else
				{
					ParticleAngle = atan2( ((*p).y()-YShift) , ((*p).x()-XShift)  );
				}

				if(theta!=0)
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
		// HEL switched on/off if turn = muliple of n
		if (SkipTurn == 0)
		{
			cout << "\n\tHEL warning: SkipTurn not set, autoset to 2" << endl;
			SkipTurn = 2;
		}
		if((Turn % SkipTurn)==0)
		{
			for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++)
			{
				//~ cout << "\n\tTurnskip Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
				if(SimpleProfile)
				{
					theta = CalcKickSimple(*p);
				}
				else
				{
					theta = CalcKickRadial(*p);
				}


				//~ if ((*p).x() < 0){	ParticleAngle = 2*pi + atan2((*p).y(), (*p).x());	}
				//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}
				//~ if ((*p).x() < 0){	ParticleAngle = pi + atan2((*p).y(), (*p).x());	}
				//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}
				if(!Elliptical)
				{
					ParticleAngle = atan2((*p).y(), (*p).x());
				}
				else
				{
					ParticleAngle = atan2( ((*p).y()-YShift) , ((*p).x()-XShift)  );
				}


				if(theta!=0)
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

}

double HollowELensProcess::GetMaxAllowedStepSize () const
{
	return currentComponent->GetLength();
}

double HollowELensProcess::CalcThetaMax (double r)
{
	if (r == 0)
	{
		return 0;
	}

	//~ cout << "L = " << EffectiveLength << " Current = " << Current << " Rigidity = " << Rigidity << " Rmax = " << Rmax << " ElectronBeta = " << ElectronBeta << endl;

	// OLD - Claiborne Smith et al. WE6RFP031 - DOESN'T WORK in the current process
	//~ ThetaMax = ((0.2 * EffectiveLength * Current) / (Rigidity * Rmax)) * ((1 + ElectronBeta)/(ElectronBeta));

	// NEW - Previtali et al. MOPWO044
	//~ ThetaMax = ( -2.0 *  EffectiveLength * Current * (1 + ElectronBeta * ProtonBeta) ) / ( 4 * pi * FreeSpacePermittivity * r * Rigidity * ElectronBeta * ProtonBeta * SpeedOfLight * SpeedOfLight);
	//~ ThetaMax = ( 2.0 *  EffectiveLength * Current * (1 + ElectronBeta * ProtonBeta) ) / ( 4 * pi * FreeSpacePermittivity * r * Rigidity * ElectronBeta * ProtonBeta * SpeedOfLight * SpeedOfLight);

	// Simplify 4 pi e0 c^2 = 1/10^-7 = 10^7
	if(ElectronDirection)
	{
		// HEL electrons travelling opposite to LHC protons (summed kick)
		ThetaMax = (2 * EffectiveLength * Current * (1 + (ElectronBeta * ProtonBeta) ) )/ ( r * 1E7 * Rigidity * ElectronBeta * ProtonBeta );
	}
	else
	{
		// HEL electrons travelling in the same direction to LHC protons (smaller kick and opposite)
		ThetaMax = -(2 * EffectiveLength * Current * (1 - (ElectronBeta * ProtonBeta) ) )/ ( r * 1E7 * Rigidity * ElectronBeta * ProtonBeta );
	}

	return ThetaMax;
}

double HollowELensProcess::CalcKickSimple (Particle &p)
{
	double thet = 0;
	double Length = 0;
	double x = 0;
	double y = 0;

	if(EffectiveLength == 0.)
	{
		cout << "HELProcess: Length = 0, setting L = 3.0[m]" << endl;
		Length = 3.0;
	}
	else
	{
		Length = EffectiveLength;
	}
	//~ cout << "\n\tHEL Length = " << Length << endl;

	// Start of HEL
	x = p.x();
	y = p.y();

	// Calculate particle transverse vector ('radius' in xy space)
	// Old for circular operation only
	//~ R = sqrt( pow((x-XOffset),2) + pow((y-YOffset),2) );
	// New takes into account elliptical operation
	R = sqrt( pow((x-XOffset+XShift),2) + pow((y-YOffset+YShift),2) );
	//~ R = sqrt( pow((x),2) + pow((y),2) );

	//~ cout << "\n\n\t\tHELProcess :: R = " << R << " Rmin = " << Rmin << endl;

	if (R < Rmin)
	{
		//~ cout << "\t\tHELProcess :: Returning zero at R = " << R << " Rmin = " << Rmin << endl;
		return 0;
	}
	else if (R < Rmax && R > Rmin)
	{
		thet = ((pow(R,2) - pow(Rmin,2))/(pow(Rmax,2) - pow(Rmin,2))) * CalcThetaMax(R);
		//thet = ((R - Rmin)/(Rmax - Rmin)) * CalcThetaMax(R);
		//cout << "\nR < Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(R) << "angle = " << ParticleAngle << endl;
		return -1* thet;

	}
	else if (R >= Rmax)
	{
		thet = CalcThetaMax(R);
		//cout << "\nR > Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(Rmax) << "angle = " << ParticleAngle << endl;
		return -1* thet;
	}
	return 0;

}

double HollowELensProcess::CalcKickSimple (double r)
{
	double thet = 0;
	double Length = 0;

	if(EffectiveLength == 0.)
	{
		cout << "HELProcess: Length = 0, setting L = 3.0[m]" << endl;
		Length = 3.0;
	}
	else
	{
		Length = EffectiveLength;
	}
	//~ cout << "\n\tHEL Length = " << Length << endl;

	// Calculate particle transverse vector ('radius' in xy space)
	R = r;

	//~ cout << "\n\n\t\tHELProcess :: R = " << R << " Rmin = " << Rmin << endl;

	if (R < Rmin)
	{
		//~ cout << "\t\tHELProcess :: Returning zero at R = " << R << " Rmin = " << Rmin << endl;
		return 0;
	}
	else if (R < Rmax && R > Rmin)
	{
		thet = ((pow(R,2) - pow(Rmin,2))/(pow(Rmax,2) - pow(Rmin,2))) * CalcThetaMax(R);
		//thet = ((R - Rmin)/(Rmax - Rmin)) * CalcThetaMax(R);
		//cout << "\nR < Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(R) << "angle = " << ParticleAngle << endl;
		return -1* thet;

	}
	else if (R >= Rmax)
	{
		thet = CalcThetaMax(R);
		//cout << "\nR > Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(Rmax) << "angle = " << ParticleAngle << endl;
		return -1* thet;
	}
	return 0;

}

double HollowELensProcess::CalcKickRadial (Particle &p)
{
	double f = 0;
	double thet = 0;
	double Length = 0;
	double x = 0;
	double y = 0;

	if(EffectiveLength == 0.)
	{
		cout << "HELProcess: Length = 0, setting L = 3.0[m]" << endl;
		Length = 3.0;
	}
	else
	{
		Length = EffectiveLength;
	}
	//~ cout << "\n\tHEL Length = " << Length << endl;

	// Start of HEL
	x = p.x();
	y = p.y();

	// Calculate particle transverse vector ('radius' in xy space)
	// Old for circular operation only
	//~ R = sqrt( pow((x-XOffset),2) + pow((y-YOffset),2) );
	// New takes into account elliptical operation
	R = sqrt( pow((x-XOffset+XShift),2) + pow((y-YOffset+YShift),2) );
	//~ R = sqrt( pow((x),2) + pow((y),2) );

	// Adapted from V. Previtali's SixTrack elense implementation

	if (R < Rmin)
	{
		return 0;
	}

	// Define boundaries between parameterisation of measured radial profile
	double r0,r1,r2,r3,r4;

	if(!LHC_Radial)
	{
		// Tevatron HEL 1.2A, 2m, 5KeV, 4-6.8sig
		r0 = (const double) 222.5;
		r1 = (const double) 252.5;
		r2 = (const double) 287;
		r3 = (const double) 364.5;
		r4 = (const double) 426.5;
	}
	else
	{
		// LHC HEL 5A, 3m, 10KeV, 4-8sig
		r1 = (const double) 265;		// 0 - initial rise				(x1.191)
		r2 = (const double) 315;		// rise - straight section		(x1.248)
		r3 = (const double) 435;		// straight - left of peak		(x1.193)
		r4 = (const double) 505;		// left - right of peak			(x1.184)
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

	double n0 = ((y1-y0)/(x1-x0)/3)*pow(x1,3)+(y0 -x0*(y1-y0)/(x1-x0))*pow(x1,2)/2 -(((y1-y0)/(x1-x0)/3)*pow(x0,3)+(y0 - x0*(y1-y0)/(x1-x0))*pow(x0,2)/2);
	double n1 = (y2-y1)/(x2-x1)*pow(x2,3)/3+(y1 - x1 *(y2-y1)/(x2-x1))*pow(x2,2)/2 -((y2-y1)/(x2-x1)*pow(x1,3)/3+(y1 - x1 *(y2-y1)/(x2-x1))*pow(x1,2)/2);
	double n2 = (y3-y2)/(x3-x2)*pow(x3,3)/3+(y2 - x2 *(y3-y2)/(x3-x2))*pow(x3,2)/2 -((y3-y2)/(x3-x2)*pow(x2,3)/3+(y2 - x2 *(y3-y2)/(x3-x2))*pow(x2,2)/2);
	double n3 = (y4-y3)/(x4-x3)*pow(x4,3)/3+(y3 - x3 *(y4-y3)/(x4-x3))*pow(x4,2)/2 -((y4-y3)/(x4-x3)*pow(x3,3)/3+(y3 - x3 *(y4-y3)/(x4-x3))*pow(x3,2)/2);
	double ntot = n0 + n1 + n2 + n3;

	if (R < x0)
	{
		f=0;
		cout << "HEL warning: Radial profile: R < x0" << endl;
	}
	else if (R < x1)
	{
		f = (((y1-y0)/(x1-x0)/3)*pow(R,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(R,2)/2-(((y1-y0)/(x1-x0)/3)*pow(x0,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(x0,2)/2))/ntot;
	}
	else if (R < x2)
	{
		f = (n0+(y2-y1)/(x2-x1)*pow(R,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(R,2)/2-((y2-y1)/(x2-x1)*pow(x1,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(x1,2)/2))/ntot;
	}
	else if (R < x3)
	{
		f = (n0+n1+(y3-y2)/(x3-x2)*pow(R,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(R,2)/2-((y3-y2)/(x3-x2)*pow(x2,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(x2,2)/2))/ntot;
	}
	else if (R < x4)
	{
		f = (n0+n1+n2+(y4-y3)/(x4-x3)*pow(R,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(R,2)/2-((y4-y3)/(x4-x3)*pow(x3,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(x3,2)/2))/ntot;
	}
	else
	{
		f = 1;
	}

	thet = -1* CalcThetaMax(R) * f;
	//cout << "\nR < Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(R) << "angle = " << ParticleAngle << endl;
	return thet;
}

double HollowELensProcess::CalcKickRadial (double r)
{
	double f = 0;
	double thet = 0;
	double Length = 0;

	if(EffectiveLength == 0.)
	{
		cout << "HELProcess: Length = 0, setting L = 3.0[m]" << endl;
		Length = 3.0;
	}
	else
	{
		Length = EffectiveLength;
	}
	//~ cout << "\n\tHEL Length = " << Length << endl;

	// Calculate particle transverse vector ('radius' in xy space)
	R = r;

	// Adapted from V. Previtali's SixTrack elense implementation

	if (R < Rmin)
	{
		return 0;
	}

	// Define boundaries between parameterisation of measured radial profile
	double r0,r1,r2,r3,r4;

	if(!LHC_Radial)
	{
		// Tevatron HEL 1.2A, 2m, 5KeV, 4-6.8sig
		r0 = (const double) 222.5;
		r1 = (const double) 252.5;
		r2 = (const double) 287;
		r3 = (const double) 364.5;
		r4 = (const double) 426.5;
	}
	else
	{
		// LHC HEL 5A, 3m, 10KeV, 4-8sig
		r1 = (const double) 265;		// 0 - initial rise				(x1.191)
		r2 = (const double) 315;		// rise - straight section		(x1.248)
		r3 = (const double) 435;		// straight - left of peak		(x1.193)
		r4 = (const double) 505;		// left - right of peak			(x1.184)
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

	double n0 = ((y1-y0)/(x1-x0)/3)*pow(x1,3)+(y0 -x0*(y1-y0)/(x1-x0))*pow(x1,2)/2 -(((y1-y0)/(x1-x0)/3)*pow(x0,3)+(y0 - x0*(y1-y0)/(x1-x0))*pow(x0,2)/2);
	double n1 = (y2-y1)/(x2-x1)*pow(x2,3)/3+(y1 - x1 *(y2-y1)/(x2-x1))*pow(x2,2)/2 -((y2-y1)/(x2-x1)*pow(x1,3)/3+(y1 - x1 *(y2-y1)/(x2-x1))*pow(x1,2)/2);
	double n2 = (y3-y2)/(x3-x2)*pow(x3,3)/3+(y2 - x2 *(y3-y2)/(x3-x2))*pow(x3,2)/2 -((y3-y2)/(x3-x2)*pow(x2,3)/3+(y2 - x2 *(y3-y2)/(x3-x2))*pow(x2,2)/2);
	double n3 = (y4-y3)/(x4-x3)*pow(x4,3)/3+(y3 - x3 *(y4-y3)/(x4-x3))*pow(x4,2)/2 -((y4-y3)/(x4-x3)*pow(x3,3)/3+(y3 - x3 *(y4-y3)/(x4-x3))*pow(x3,2)/2);
	double ntot = n0 + n1 + n2 + n3;

	if (R < x0)
	{
		f=0;
		cout << "HEL warning: Radial profile: R < x0" << endl;
	}
	else if (R < x1)
	{
		f = (((y1-y0)/(x1-x0)/3)*pow(R,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(R,2)/2-(((y1-y0)/(x1-x0)/3)*pow(x0,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(x0,2)/2))/ntot;
	}
	else if (R < x2)
	{
		f = (n0+(y2-y1)/(x2-x1)*pow(R,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(R,2)/2-((y2-y1)/(x2-x1)*pow(x1,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(x1,2)/2))/ntot;
	}
	else if (R < x3)
	{
		f = (n0+n1+(y3-y2)/(x3-x2)*pow(R,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(R,2)/2-((y3-y2)/(x3-x2)*pow(x2,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(x2,2)/2))/ntot;
	}
	else if (R < x4)
	{
		f = (n0+n1+n2+(y4-y3)/(x4-x3)*pow(R,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(R,2)/2-((y4-y3)/(x4-x3)*pow(x3,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(x3,2)/2))/ntot;
	}
	else
	{
		f = 1;
	}

	thet = -1* CalcThetaMax(R) * f;
	//cout << "\nR < Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(R) << "angle = " << ParticleAngle << endl;
	return thet;
}

void HollowELensProcess::SetRadii (double rmin, double rmax)
{
	g = rmax/rmin;
	cout << "\n\tHEL warning: HEL radii not set using beam envelope, and not aligned to beam orbit" << endl;
	Rmin = rmin;
	Rmax = rmax;
	Rmin_original = Rmin;
	Rmax_original = Rmax;
}

void HollowELensProcess::SetRadiiSigma (double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
{
	g = rmax / rmin;
	//How many HELs in lattice?
	int Hel_no = 0;
	//Element no of last HEL
	int Hel_ID = 0;

	bool find_HEL_no = 1;
	if (find_HEL_no)
	{
		vector<AcceleratorComponent*> Elements;
		size_t n_collimators = model->ExtractTypedElements(Elements,"*");
		size_t nelm = Elements.size();

		for (size_t n = 0; n < nelm; n++)
		{
			if(Elements[n]->GetType() == "HollowElectronLens")
			{
				Hel_no++;
				Hel_ID = n;
			}
		}

		if(Hel_no == 0)
		{
			cout << "Error: No HEL found in HollowELensProcess::SetRadiiSigma " << endl;
		}
		else if(Hel_no == 1)
		{
			cout << "HollowELensProcess::SetRadiiSigma : 1 HEL found at element No " << Hel_ID <<  endl;
		}
		else if (Hel_no > 1)
		{
			cout << "HollowELensProcess::SetRadiiSigma : More than 1 HELs found, last at element No " << Hel_ID <<  endl;
		}
	}

	vector<HollowElectronLens*> HELs;
	size_t n_Hels = model->ExtractTypedElements(HELs,"*");

	cout << " HollowELensProcess::SetRadiiSigma : find_HEL : found " << HELs.size() << " HELs" << endl;
	//Calculate sigma to adjust radii as required

	double sigma_x = 0;
	double beta_x = 0;
	double alpha_x = 0;
	double gamma_x = 0;
	double sigma_xp = 0;

	double sigma_y = 0;
	double beta_y = 0;
	double alpha_y = 0;
	double gamma_y = 0;
	double sigma_yp = 0;

	for (vector<HollowElectronLens*>::iterator it = HELs.begin(); it!= HELs.end(); it++)
	{

		for(int j = 0; j <= twiss->NumberOfRows(); j++)
		{
			if((*it)->GetComponentLatticePosition() == twiss->Value(0,0,0,j))
			{
				// Note that for a thin lens (zero length) HEL, this will be true twice
				// as the element prior to the drift has the same S.
				// As the HEL will always be the second element with the same S
				// this doesn't cause any issues. If however there is a zero length element after
				// the HEL, the values below should be equal as well.

				cout << " S_HEL = " << (*it)->GetComponentLatticePosition() << "m" << endl;
				cout << " S_twiss = " << twiss->Value(0,0,0,j) << "m" << endl;

				cout << "HollowELensProcess::SetRadiiSigma : j value = " << j << endl;

				//Note that this is currently only a horizontal HEL
				beta_x = twiss->Value(1,1,1,j);		//Beta x
				beta_y = twiss->Value(3,3,2,j);		//Beta y

				XOffset = twiss->Value(1,0,0,j);
				YOffset = twiss->Value(3,0,0,j);
				twiss->PrintTable(std::cout,j-1,j);
				twiss->PrintTable(std::cout,j+1,j+2);

				//~ double sigma_y = sqrt(beta_y * emittance_y);
				sigma_x = sqrt(beta_x * emittance_x);
				sigma_y = sqrt(beta_y * emittance_y);

				alpha_x = -1*twiss->Value(1,2,1,j);
				gamma_x = (1 + (alpha_x*alpha_x) )/ beta_x;

				alpha_y = -1*twiss->Value(3,4,2,j);
				gamma_y = (1 + (alpha_y*alpha_y) )/ beta_y;

				sigma_xp = sqrt(gamma_x * emittance_x);
				sigma_yp = sqrt(gamma_y * emittance_y);

			}
		}
	}

	cout << "HollowELensProcess::SetRadiiSigma : Beta_x = " << beta_x << " Sigma_x = " << sigma_x << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_x = " << -alpha_x << " Sigma_xp = " << sigma_xp << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Beta_y = " << beta_y << " Sigma_y = " << sigma_y << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_y = " << -alpha_y << " Sigma_yp = " << sigma_yp << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Offset_x = " << XOffset << " Offset_y = " << YOffset << endl;

	Rmin = rmin * sigma_x;
	Rmax = rmax * sigma_x;
	Rmin_original = Rmin;
	Rmax_original = Rmax;
	Sigma_x = sigma_x;
	Sigma_y = sigma_y;

	cout << "HollowELensProcess::SetRadiiSigma : RMax = " << Rmax << " RMin= " << Rmin << endl;

}

void HollowELensProcess::SetRadiiSigma (double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss, double P0)
{

	g = rmax/rmin;
	//How many HELs in lattice?
	int Hel_no = 0;
	//Element no of last HEL
	int Hel_ID = 0;

	bool find_HEL_no = 1;
	if (find_HEL_no)
	{
		vector<AcceleratorComponent*> Elements;
		size_t n_collimators = model->ExtractTypedElements(Elements,"*");
		size_t nelm = Elements.size();

		for (size_t n = 0; n < nelm; n++)
		{
			if(Elements[n]->GetType() == "HollowElectronLens")
			{
				Hel_no++;
				Hel_ID = n;
			}
		}

		if(Hel_no == 0)
		{
			cout << "Error: No HEL found in HollowELensProcess::SetRadiiSigma " << endl;
		}
		else if(Hel_no == 1)
		{
			cout << "HollowELensProcess::SetRadiiSigma : 1 HEL found at element No " << Hel_ID <<  endl;
		}
		else if (Hel_no > 1)
		{
			cout << "HollowELensProcess::SetRadiiSigma : More than 1 HELs found, last at element No " << Hel_ID <<  endl;
		}
	}

	vector<HollowElectronLens*> HELs;
	size_t n_Hels = model->ExtractTypedElements(HELs,"*");

	cout << " HollowELensProcess::SetRadiiSigma : find_HEL : found " << HELs.size() << " HELs" << endl;
	//Calculate sigma to adjust radii as required

	double sigma_x = 0;
	double beta_x = 0;
	double alpha_x = 0;
	double gamma_x = 0;
	double sigma_xp = 0;

	double sigma_y = 0;
	double beta_y = 0;
	double alpha_y = 0;
	double gamma_y = 0;
	double sigma_yp = 0;

	int stored_j = 0;
	double stored_s = 0;

	for (vector<HollowElectronLens*>::iterator it = HELs.begin(); it!= HELs.end(); it++)
	{

		for(int j = 0; j <= twiss->NumberOfRows(); j++)
		{
			if((*it)->GetComponentLatticePosition() == twiss->Value(0,0,0,j))
			{
				// Note that for a thin lens (zero length) HEL, this will be true twice
				// as the element prior to the drift has the same S.
				// As the HEL will always be the second element with the same S
				// this doesn't cause any issues. If however there is a zero length element after
				// the HEL, the values below should be equal as well.

				stored_j = j;
				stored_s = (*it)->GetComponentLatticePosition();
				cout << " S_HEL = " << (*it)->GetComponentLatticePosition() << "m" << endl;
				cout << " S_twiss = " << twiss->Value(0,0,0,j) << "m" << endl;

				cout << "HollowELensProcess::SetRadiiSigma : j value = " << j << endl;

				//Note that this is currently only a horizontal HEL
				beta_x = twiss->Value(1,1,1,j);		//Beta x
				beta_y = twiss->Value(3,3,2,j);		//Beta y

				XOffset = twiss->Value(1,0,0,j);
				YOffset = twiss->Value(3,0,0,j);
				twiss->PrintTable(std::cout,j-1,j);
				twiss->PrintTable(std::cout,j+1,j+2);

				//~ double sigma_y = sqrt(beta_y * emittance_y);
				sigma_x = sqrt(beta_x * emittance_x);
				sigma_y = sqrt(beta_y * emittance_y);

				alpha_x = -1*twiss->Value(1,2,1,j);
				gamma_x = (1 + (alpha_x*alpha_x) )/ beta_x;

				alpha_y = -1*twiss->Value(3,4,2,j);
				gamma_y = (1 + (alpha_y*alpha_y) )/ beta_y;

				sigma_xp = sqrt(gamma_x * emittance_x);
				sigma_yp = sqrt(gamma_y * emittance_y);

			}
		}
	}

	PhaseAdvance* PA = new PhaseAdvance(model, twiss, P0);
	double Mux = PA->GetPhaseAdvanceX(Hel_ID);
	double Muy = PA->GetPhaseAdvanceY(Hel_ID);

	cout << "\nHollowELensProcess::SetRadiiSigma : S = " << stored_s << endl;
	cout << "\nHollowELensProcess::SetRadiiSigma : Beta_x = " << beta_x << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Beta_y = " << beta_y << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_x = " << -alpha_x << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_y = " << -alpha_y << endl;
	cout << "\nHollowELensProcess::SetRadiiSigma : Sigma_x = " << sigma_x << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Sigma_y = " << sigma_y << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Sigma_xp = " << sigma_xp << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Sigma_yp = " << sigma_yp << endl;
	cout << "\nHollowELensProcess::SetRadiiSigma : Offset_x = " << XOffset << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Offset_y = " << YOffset << endl;
	cout << "\nHollowELensProcess::SetRadiiSigma : mu_x = " << Mux << endl;
	cout << "HollowELensProcess::SetRadiiSigma : mu_y = " << Muy << endl;

	Rmin = rmin * sigma_x;
	Rmax = rmax * sigma_x;
	Rmin_original = Rmin;
	Rmax_original = Rmax;
	Sigma_x = sigma_x;
	Sigma_y = sigma_y;

	cout << "HollowELensProcess::SetRadiiSigma : RMax = " << Rmax << " RMin= " << Rmin << endl;

}

void HollowELensProcess::OutputProfile(std::ostream* os, double E, double min, double max)
{
	// Have to set ProtonBeta if running before it is set
	double Gamma_p = LorentzGamma(E, ProtonMass);
	ProtonBeta = LorentzBeta(Gamma_p);

	cout << " Rmin = " << Rmin << ", = " << Rmin/Sigma_x << " sigma" << endl;
	cout << " Rmax = " << Rmax << ", = " << Rmax/Sigma_x << " sigma" << endl;
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
		(*os) << (r/Sigma_x) <<"\t"<< CalcKickRadial(r) <<"\t"<< CalcKickSimple(r) <<"\t"<< sqrt(pow(CalcKickRadial(r),2)) <<"\t"<< sqrt(pow(CalcKickSimple(r),2)) << endl;
		//~ (*os) << (r/Sigma_x) <<"\t"<< CalcKickRadial(r) <<"\t"<< CalcKickSimple(r)  << endl;
		//~ (*os) << (r/Sigma_x) <<"\t"<< sqrt(pow(CalcKickRadial(r),2)) <<"\t"<< sqrt(pow(CalcKickSimple(r),2)) << endl;
		r += ((max-min)/points) * Sigma_x;
	}
}

void HollowELensProcess::OutputFootprint(std::ostream* os, int npart)
{
	cout << " Rmin = " << Rmin << ", = " << Rmin/Sigma_x << " sigma" << endl;
	cout << " Rmax = " << Rmax << ", = " << Rmax/Sigma_x << " sigma" << endl;
	cout << " L = " << EffectiveLength << endl;
	cout << " Current = " << Current << endl;
	cout << " Brho = " << Rigidity << endl;
	cout << " ElectronBeta = " << ElectronBeta << endl;
	cout << " ProtonBeta = " << ProtonBeta << endl;

	// Make a 'distribution' in xy phase space
	double x_probe;
	double y_probe;

	// using up to 3*R
	double max_xy = 3*(sqrt(pow(Rmax,2))/2);

	(*os) << "#r\tx\ty\tr" << endl;

	// Note here we subtract XShift, for the particle distribution we add XShift
	for(int i = 0; i<npart; ++i)
	{
		x_probe = RandomNG::uniform(-max_xy,max_xy) - XOffset - XShift;
		y_probe = RandomNG::uniform(-max_xy,max_xy) - YOffset - YShift;
		if( ((sqrt( pow((x_probe-XOffset-XShift),2) + pow((y_probe-YOffset-YShift),2) )) > Rmin) &&((sqrt( pow((x_probe-XOffset-XShift),2) + pow((y_probe-YOffset-YShift),2) )) < Rmax) )
		{
			(*os) << x_probe << "\t" << y_probe << "\t" << sqrt(pow(x_probe,2) + pow(y_probe,2)) << endl;
		}
	}
}

void HollowELensProcess::SetElectronDirection(bool dir)
{
	ElectronDirection = dir;
	if(ElectronDirection)
	{
		cout << "HELProcess: electrons travelling opposite to protons: negative (focussing) kick" << endl;
	}
	else
	{
		cout << "HELProcess: electrons travelling in the same direction as protons: positive (defocussing) kick" << endl;
	}
}

void HollowELensProcess::SetEllipticalMatching(bool io)
{

	Elliptical = io;

	// Radii in sigma
	double rmin = Rmin_original/Sigma_x;
	double rmax = Rmax_original/Sigma_x;

	cout << "\n\t HELProcess : Using Elliptical Operation" << endl;

	// Need to set SemiMajor and SemiMinor axes for our ellipse - will be used to fnd HEL R_min
	// The ratio of Rmin to Max is set by hardware and is 2 for the LHC, and shoulf be set in SetRadiiSigma or SetRadii
	if(g ==0.)
	{
		g = 2.0;
	}
	// Also need to set XShift or YShift for co-ordinate transforms
	if(Elliptical)
	{
		// x larger than y, shift co-ordinates up in y
		if(Sigma_x > Sigma_y)
		{
			SemiMinor = Sigma_y * rmin;
			SemiMajor = Sigma_x * rmin;
			Rmin = sqrt(SemiMajor/SemiMinor)*( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
			//~ Rmin = ( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
			Rmax = g*Rmin;
			YShift = SemiMinor - Rmin;
		}
		// y larger than x, shift co-ordinates right in x (i.e. assuming beam1 for LHC)
		else if(Sigma_x < Sigma_y)
		{
			SemiMinor = Sigma_x * rmin;
			SemiMajor = Sigma_y * rmin;
			Rmin = sqrt(SemiMajor/SemiMinor)*( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
			//~ Rmin = ( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
			Rmax = g*Rmin;
			XShift = SemiMinor - Rmin;
		}

		EllipticalSet = 1;
	}
	else
	{
		EllipticalSet = 0;
	}
}

void HollowELensProcess::EllipticalAdjust()
{

	// Radii in sigma
	double rmin = Rmin_original/Sigma_x;
	double rmax = Rmax_original/Sigma_x;

	// x larger than y, shift co-ordinates up in y
	if(Sigma_x > Sigma_y)
	{
		SemiMinor = Sigma_y * rmin;
		SemiMajor = Sigma_x * rmin;
		Rmin = sqrt(SemiMajor/SemiMinor)*( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
		//~ Rmin = ( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
		Rmax = g*Rmin;
		YShift = SemiMinor - Rmin;
	}
	// y larger than x, shift co-ordinates right in x (i.e. assuming beam1 for LHC)
	else if(Sigma_x < Sigma_y)
	{
		SemiMinor = Sigma_x * rmin;
		SemiMajor = Sigma_y * rmin;
		Rmin = sqrt(SemiMajor/SemiMinor)*( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
		//~ Rmin = ( pow(SemiMajor,2) + pow(SemiMinor,2) ) / (2 * SemiMinor);
		Rmax = g*Rmin;
		XShift = SemiMinor - Rmin;
	}

}

void HollowELensProcess::EllipticalAdjust(int compass)
{

	// Radii in sigma
	double rmin = Rmin_original/Sigma_x;
	double rmax = Rmax_original/Sigma_x;
	bool horizontal = 0;

	// Radii in Sigma
	if(Sigma_x > Sigma_y)
	{
		SemiMinor = Sigma_y * rmin;
		SemiMajor = Sigma_x * rmin;
		horizontal = 0;
	}
	else
	{
		SemiMinor = Sigma_x * rmin;
		SemiMajor = Sigma_y * rmin;
		horizontal = 1;
	}
}

}
