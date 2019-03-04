/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef HollowElectronLens_h
#define HollowElectronLens_h 1

#include "merlin_config.h"
#include "Drift.h"

class ComponentTracker;

/**
 * HEL operation modes
 */
typedef enum
{
	DC,
	AC,
	Diffusive,
	Turnskip

} OperationMode;

/**
 * A HollowElectronLens provides a method of active halo control.
 * It can be used as a soft scraper. Optically the element is a drift.
 * See HollowELensProcess.
 */

class HollowElectronLens: public Drift
{
public:

	HollowElectronLens(const string& id, double len, int mode, double current, double beta_e, double rigidity, double
		length_e);

	/**
	 *	Returns the type string for this component.
	 *	@return Component type string
	 */
	virtual const string& GetType() const;

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 *	@return Accelerator component class unique index
	 */
	virtual int GetIndex() const;

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Rotates the component 180 degrees about its local Y axis.
	 */
	virtual void RotateY180();

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	virtual double GetRmax() const
	{
		return Rmax;
	}

	virtual double GetRmin() const
	{
		return Rmin;
	}

	virtual void SetRmax(double rmax);
	virtual void SetRmin(double rmin);
	void SetRadii(double rmin, double rmax);

	/**
	 * Set the effective length of the e- lens
	 */
	virtual void SetEffectiveLength(double l_e)
	{
		EffectiveLength = l_e;
	}

	/**
	 * Set the type of HEL operation required
	 */
	virtual void SetOpMode(OperationMode mode)
	{
		OMode = mode;
	}

	/**
	 * Set variables for AC mode operation
	 */
	virtual void SetAC(double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi);

	/**
	 * SetTurnSkip
	 */
	virtual void SetTurnskip(int skip);

	/**
	 * Change to radial (measured) profile, simple (perfect) is default
	 */
	virtual void SetRadialProfile()
	{
		SimpleProfile = 0;
	}
	virtual void SetPerfectProfile()
	{
		SimpleProfile = 1;
	}
	virtual void SetLHCRadialProfile()
	{
		LHC_Radial = 1;
	}

	/**
	 * Change electron direction (default opposite protons = 1)
	 */
	virtual void SetElectronDirection(bool dir);

	/**
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;

private:
	double Rmin;
	double Rmax;

public:
	// Hardware parameters
	double Current;
	double ElectronBeta;
	double Rigidity;

	double EffectiveLength;
	double Sigma_x;

	// Variables
	double XOffset;
	double YOffset;

	//For AC mode
	double Tune;
	double DeltaTune;
	double TuneVarPerStep;
	double TurnsPerStep;
	double Multiplier;
	double Nstep;

	int Turn;
	int SkipTurn;

	double MinTune;

	bool ACSet;                 /// AC mode variable set?
	bool SimpleProfile;         /// 1 = use perfect HEL profile, 0 = use parameterisation of measured LHC prototype cathode profile
	bool AlignedToOrbit;        /// Is the HEL aligned to the closed orbit
	bool ElectronDirection;     /// 1 = opposite protons (-ve kick), 0 = same as protons (+ve kick)
	bool LHC_Radial;            /// 1 = use empirically scaled measured radial profile (LHC hardware), 0 = use measured radial HEL profile

	OperationMode OMode;
};

#endif
