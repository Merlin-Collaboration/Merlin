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
// Last Edited: 19.09.15 HR
// 
/////////////////////////////////////////////////////////////////////////
#include <cmath>

#include "NumericalUtils/utils.h"

#include "Random/RandomNG.h"

#include "BeamDynamics/ParticleTracking/HollowELensProcess.h"

#include "RingDynamics/LatticeFunctions.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace std;

namespace ParticleTracking {


HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity), ACSet(0), SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0)
{
	if (mode == 0){OMode = DC;}
	else if (mode == 1){OMode = AC;}
	else if (mode == 2){OMode = Diffusive;}
	else if (mode == 3){OMode = Turnskip;}
	else { cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;}
}

HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity), EffectiveLength(length_e), ACSet(0), SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0)
{
	if (mode == 0){OMode = DC;}
	else if (mode == 1){OMode = AC;}
	else if (mode == 2){OMode = Diffusive;}
	else if (mode == 3){OMode = Turnskip;}
	else { cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;}
	
}

HollowELensProcess::HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
	: ParticleBunchProcess("HOLLOW ELECTRON LENS", priority), Current(current), ElectronBeta(beta_e), Rigidity(rigidity), ACSet(0), SimpleProfile(1), AlignedToOrbit(0), XOffset(0), YOffset(0), Turn(0), SkipTurn(0)
{
	if (mode == 0){OMode = DC;}
	else if (mode == 1){OMode = AC;}
	else if (mode == 2){OMode = Diffusive;}
	else if (mode == 3){OMode = Turnskip;}
	else { cout << "\tHEL operation mode invalid. Please choose between: \n\t int 0 = DC \n\t int 1 = AC \n\t int 2 = Diffusive \n\t int 3 = Turnskip" << endl;}
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
	active = (currentBunch!=0) && (aHollowELens);
	
	if(active) {
	    currentComponent = &component;
	    EffectiveLength = currentComponent->GetLength();		
		
		double Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
		}
		else {
		    currentComponent = 0;
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
	double theta = 0;
	double Gamma_p = 0;
	//cout << "Entered Hollow Electron Lens Process" << endl;

	ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(), currentBunch->GetTotalCharge()/currentBunch->size());
	newbunch->clear();
	newbunch->swap(*currentBunch);

	if(ProtonBeta == 0){
		Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		ProtonBeta = LorentzBeta(Gamma_p);
	}	
	
	// Have to increment Turn as the process doesn't have access to the turn value from user code
	++Turn;

	switch (OMode){
		case DC:{
			//HEL always on
			for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++) 
			{				
				//~ cout << "\n\tDC Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;  
				if(SimpleProfile){theta = CalcKickSimple(*p);}
				else{theta = CalcKickRadial(*p);}
							
				ParticleAngle = atan2((*p).y(), (*p).x());
				//~ cout << "\n\tX = " << (*p).x() << endl;
				//~ cout << "\n\tY = " << (*p).y() << endl;
				//~ cout << "\n\tAngle = " << ParticleAngle << endl;
				//~ cout << "\n\tRadius = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2)) << endl;
				//~ cout << "\n\tRadius = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2))/ sqrt(pow(293.031E-6,2) + pow(287.801E-6,2))<< " sigma" << endl;
				
				if(theta!=0){
        
					//~ cout << " x = " << (*p).x() << endl;
					//~ cout << " y = " << (*p).y() << endl;
					//~ cout << " xp = " << (*p).xp() << endl;
					//~ cout << " yp = " << (*p).yp() << endl;
					//~ cout << " ParticleAngle = " << ParticleAngle << endl;
					//~ cout << " R = " << sqrt(pow((*p).x(),2) + pow((*p).y(),2)) << endl;
					//~ cout << " Rmin = " << Rmin << endl;
					//~ cout << " Rmax = " << Rmax << endl;
					//~ cout << " L = " << EffectiveLength << endl;
					//~ cout << " max_kick = " << CalcThetaMax(R) << endl;
					//~ cout << " Current = " << Current << endl;
					//~ cout << " Brho = " << Rigidity << endl;
					//~ cout << " Gamma_p = " << Gamma_p << endl;
					//~ cout << " theta = " << theta << endl;
					//~ cout << "\n" << endl;
					
					
					//~ // Particle phase space angle and amplitude (radius)
					(*p).xp() -= theta * cos(ParticleAngle);			
					(*p).yp() -= theta * sin(ParticleAngle);
				}
			}
		}
		break;
		case AC:{
			//~ cout << "\n\tAC Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
			// Resonant HEL kick - Adapted from V. Previtali's SixTrack elense
			if(ACSet){
				for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++) 
				{
					if(SimpleProfile){theta = CalcKickSimple(*p);}
					else{theta = CalcKickRadial(*p);}
					
					if(theta!=0){
						
						if( (TuneVarPerStep !=0) && (DeltaTune !=0) ){						
							OpTune = MinTune + fmod((floor(Turn/TurnsPerStep)),(Nstep));
						}
						else{
							OpTune = Tune;
						}
						
						Phi = Multiplier * ( Turn * OpTune * 2 * pi );
						theta *= 0.5*(1 + cos(Phi));
										
						ParticleAngle = atan2((*p).y(), (*p).x());
						
						// Particle phase space angle and amplitude (radius)
						(*p).xp() -= theta * cos(ParticleAngle);			
						(*p).yp() -= theta * sin(ParticleAngle);
					}	
				}
			}
			else{cout << "\n\tHEL Warning: AC variables not set" << endl;}
		}	
		break;
		case Diffusive:
		{
			//~ cout << "\n\tDiffusive Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
			// HEL randomly switched on/off on a turn by turn basis
			double rando = RandomNG::uniform(-1,1);
			
			if (rando >=0){
				for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++) 
				{
					if(SimpleProfile){theta = CalcKickSimple(*p);}
					else{theta = CalcKickRadial(*p);}
												
				
					//~ if ((*p).x() < 0){	ParticleAngle = pi + atan2((*p).y(), (*p).x());	}
					//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}		
					//~ if ((*p).x() < 0){	ParticleAngle = 2*pi + atan2((*p).y(), (*p).x());	}
					//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}	
					ParticleAngle = atan2((*p).y(), (*p).x());	
						
					if(theta!=0){
						// Particle phase space angle and amplitude (radius)
						(*p).xp() -= theta * cos(ParticleAngle);			
						(*p).yp() -= theta * sin(ParticleAngle);
					}			

				}
			}
		}
		break;
		case Turnskip:
		{
			// HEL switched on/off if turn = muliple of n
			if (SkipTurn == 0){cout << "\n\tHEL warning: SkipTurn not set, autoset to 2" << endl; SkipTurn = 2;}
			if((Turn % SkipTurn)==0){
				for(PSvectorArray::iterator p = newbunch->begin(); p!=newbunch->end(); p++) 
				{
					//~ cout << "\n\tTurnskip Mode, Rmin = " << Rmin << " Rmax = " << Rmax << endl;
					if(SimpleProfile){theta = CalcKickSimple(*p);}
					else{theta = CalcKickRadial(*p);}
												
			
					//~ if ((*p).x() < 0){	ParticleAngle = 2*pi + atan2((*p).y(), (*p).x());	}
					//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}		
					//~ if ((*p).x() < 0){	ParticleAngle = pi + atan2((*p).y(), (*p).x());	}
					//~ else{				ParticleAngle = atan2((*p).y(), (*p).x());		}
					ParticleAngle = atan2((*p).y(), (*p).x());		
	
					
					if(theta!=0){
						// Particle phase space angle and amplitude (radius)
						(*p).xp() -= theta * cos(ParticleAngle);			
						(*p).yp() -= theta * sin(ParticleAngle);
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
	if (r == 0){return 0;}

	//cout << "L = " << EffectiveLength << " Current = " << Current << " Rigidity = " << Rigidity << " Rmax = " << Rmax << " ElectronBeta = " << ElectronBeta << endl;

	//OLD - Claiborne Smith et al. WE6RFP031
	//ThetaMax = ((0.2 * EffectiveLength * Current) / (Rigidity * Rmax)) * ((1 + ElectronBeta)/(ElectronBeta));
	
	//NEW - Previtali et al. MOPWO044
	//ThetaMax = ( -2.0 *  EffectiveLength * Current * (1 + ElectronBeta * ProtonBeta) ) / ( 4 * pi * FreeSpacePermittivity * r * Rigidity * ElectronBeta * ProtonBeta * SpeedOfLight * SpeedOfLight);  

	//Simplify 4 pi e0 c^2 = 1/10^-7 = 10^7	
	ThetaMax = (2 * EffectiveLength * Current * (1 + (ElectronBeta * ProtonBeta) ) )/ ( r * 1E7 * Rigidity * ElectronBeta * ProtonBeta );	

	return ThetaMax;
}

double HollowELensProcess::CalcKickSimple (Particle &p)
{
	double thet = 0;
	double Length = 0;
	double x = 0;
	double y = 0;
	
	if(EffectiveLength == 0.){cout << "HELProcess: Length = 0, setting L = 2.0[m]" << endl; Length = 2.0;}
	else{Length = EffectiveLength;}
	//~ cout << "\n\tHEL Length = " << Length << endl;	

	// Start of HEL
	x = p.x();
	y = p.y();	
	
	// Calculate particle transverse vector ('radius' in xy space)
	R = sqrt( pow((x-XOffset),2) + pow((y-YOffset),2) );
	//~ R = sqrt( pow((x),2) + pow((y),2) );
	
	//~ cout << "\n\n\t\tHELProcess :: R = " << R << " Rmin = " << Rmin << endl; 

	if (R < Rmin) {
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
	
	if(EffectiveLength == 0.){Length = 2.0;}
	else{Length = EffectiveLength;}
	//~ cout << "\n\tHEL Length = " << Length << endl;	

	// Start of HEL
	x = p.x();
	y = p.y();	
	
	// Calculate particle transverse vector ('radius' in xy space)
	R = sqrt( pow((x-XOffset),2) + pow((y-YOffset),2) );
	//~ R = sqrt( pow((x),2) + pow((y),2) );
	
	// Adapted from V. Previtali's SixTrack elense implementation 
	
	if (R < Rmin) {
		return 0;
	}
	
	// Define boundaries between parameterisation of measured radial profile
	const double r0 = 222.5;
	const double r1 = 252.5;
	const double r2 = 287;
	const double r3 = 364.5;
	const double r4 = 426.5;
	
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
    
    if (R < x0){
                 f=0;
                 cout << "HEL warning: Radial profile: R < x0" << endl;     
				}
	else if (R < x1){
                 f = (((y1-y0)/(x1-x0)/3)*pow(R,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(R,2)/2-(((y1-y0)/(x1-x0)/3)*pow(x0,3)+(y0-x0 *(y1-y0)/(x1-x0))*pow(x0,2)/2))/ntot;
				}
    else if (R < x2){
                 f = (n0+(y2-y1)/(x2-x1)*pow(R,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(R,2)/2-((y2-y1)/(x2-x1)*pow(x1,3)/3+(y1-x1*(y2-y1)/(x2-x1))*pow(x1,2)/2))/ntot;
	}
    else if (R < x3){
                f = (n0+n1+(y3-y2)/(x3-x2)*pow(R,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(R,2)/2-((y3-y2)/(x3-x2)*pow(x2,3)/3+(y2-x2*(y3-y2)/(x3-x2))*pow(x2,2)/2))/ntot;
	}
    else if (R < x4){
                f = (n0+n1+n2+(y4-y3)/(x4-x3)*pow(R,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(R,2)/2-((y4-y3)/(x4-x3)*pow(x3,3)/3+(y3-x3*(y4-y3)/(x4-x3))*pow(x3,2)/2))/ntot;
	}
    else{ 
                f = 1;                 
    } 
	
	thet = -1* CalcThetaMax(R) * f;
	//cout << "\nR < Rmax, theta = " << thet << " thetamax = " << CalcThetaMax(R) << "angle = " << ParticleAngle << endl;
	return thet;
}


void HollowELensProcess::SetRadii (double rmin, double rmax)
{
	cout << "\n\tHEL warning: HEL radii not set using beam envelope, and not aligned to beam orbit" << endl;
	Rmin = rmin;
	Rmax = rmax;
}

void HollowELensProcess::SetRadiiSigma (double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
{
	
	//How many HELs in lattice?
	int Hel_no = 0;
	//Element no of last HEL
	int Hel_ID = 0;
	
	bool find_HEL_no = 1;
	if (find_HEL_no){
		vector<AcceleratorComponent*> Elements;
		size_t n_collimators = model->ExtractTypedElements(Elements,"*");
		size_t nelm = Elements.size();
		
		for (size_t n = 0; n < nelm; n++){
			if(Elements[n]->GetType() == "HollowElectronLens"){			
				Hel_no++;
				Hel_ID = n;			
			}
		}
		
		if(Hel_no == 0){ cout << "Error: No HEL found in HollowELensProcess::SetRadiiSigma " << endl;}
		else if(Hel_no == 1){ cout << "HollowELensProcess::SetRadiiSigma : 1 HEL found at element No " << Hel_ID <<  endl; }
		else if (Hel_no > 1){ cout << "HollowELensProcess::SetRadiiSigma : More than 1 HELs found, last at element No " << Hel_ID <<  endl;}
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
	
	for (vector<HollowElectronLens*>::iterator it = HELs.begin(); it!= HELs.end(); it++){
		
		for(int j = 0; j <= twiss->NumberOfRows(); j++){
			if((*it)->GetComponentLatticePosition() == twiss->Value(0,0,0,j)){
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
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_x = " << alpha_x << " Sigma_xp = " << sigma_xp << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Beta_y = " << beta_y << " Sigma_y = " << sigma_y << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Alpha_y = " << alpha_y << " Sigma_yp = " << sigma_yp << endl;
	cout << "HollowELensProcess::SetRadiiSigma : Offset_x = " << XOffset << " Offset_y = " << YOffset << endl;
	
	Rmin = rmin * sigma_x;
	Rmax = rmax * sigma_x;
	
	cout << "HollowELensProcess::SetRadiiSigma : RMax = " << Rmax << " RMin= " << Rmin << endl;
	
}

double HollowELensProcess::CalcSigma (double emittance_x, double emittance_y, LatticeFunctionTable* twiss)
{
	//TODO necessary?
	return 0;
}

}; // end namespace ParticleTracking
