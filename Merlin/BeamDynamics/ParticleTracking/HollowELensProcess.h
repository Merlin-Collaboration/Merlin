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
#ifndef HollowELensProcess_h
#define HollowELensProcess_h 1

#include "merlin_config.h"

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/StdComponent/HollowElectronLens.h"

#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

#include "RingDynamics/LatticeFunctions.h"

namespace ParticleTracking {

// HEL operation modes
typedef enum {DC, AC, Diffusive, Turnskip} OperationMode;


class HollowELensProcess : public ParticleBunchProcess
{
public:
    //	Constructor 
    HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity);
    HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e);
    HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss);

    //	Initialise this process with the specified Bunch. If
    //	bunch is not a ParticleBunch object, the process becomes
    //	inactive.
    virtual void InitialiseProcess (Bunch& bunch);

    //	Sets the current accelerator component.
    virtual void SetCurrentComponent (AcceleratorComponent& component);

    //	Preform the process for the specified step ds.
    virtual void DoProcess (double ds);

    //	Returns the current maximum step length for this process.
    virtual double GetMaxAllowedStepSize () const;
    
 	// Set minimum and maximum e- beam radii in [m] or [sigma]
    virtual void SetRadii (double rmin, double rmax);
    virtual void SetRadiiSigma (double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss);
    
    // Used to align the e- lens with the beam envelope
    virtual double CalcSigma (double emittance_x, double emittance_y, LatticeFunctionTable* twiss);
    
    // Set the effective length of the e- lens
    virtual void SetEffectiveLength (double l_e) {EffectiveLength = l_e;}
    
    // Calculates the theta kick given by the e- lens
    virtual double CalcThetaMax (double r);

	// Set the type of HEL operation required
	virtual void SetOpMode (OperationMode mode){OMode = mode;}

	// Set variables for AC mode operation
	virtual void SetAC (double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi);
	
	//SetTurnSkip
	virtual void SetTurnskip (int skip);
	
	// Use simple profile to calculate kick
	virtual double CalcKickSimple (Particle &p);
	
	// Use radial (measured) profile to calculate kick
	virtual double CalcKickRadial (Particle &p);

	// Change to radial (measured) profile, simple (perfect) is default
	virtual void SetRadialProfile(){SimpleProfile = 0;}
	virtual void SetPerfectProfile(){SimpleProfile = 1;}
	
	virtual void OutputKick(std::ostream* os){}


private:
    // Data Members for Class Attributes
	double Current;
	double ElectronBeta;
	double Rigidity;
	double ProtonBeta;
	double EffectiveLength;
	
	double ThetaMax;
	double ParticleAngle;
	double R;
	double Rmin;
	double Rmax;

	bool ACSet;
	bool SimpleProfile;
	bool AlignedToOrbit;
	double XOffset;
	double YOffset;

	//For AC mode
	double Tune;
	double DeltaTune;
	double TuneVarPerStep;
	double TurnsPerStep;
	double Multiplier;	
	double Nstep;
	double OpTune;
	double Phi;
	int Turn;
	int SkipTurn;

	double MinTune;
	double MaxTune;

	OperationMode OMode;

};


} // end namespace ParticleTracking
#endif
