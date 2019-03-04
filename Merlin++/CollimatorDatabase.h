/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _COLLIMATOR_DATABASE_H_
#define _COLLIMATOR_DATABASE_H_

#include <fstream>
#include <string>
#include <vector>

#include "AcceleratorModel.h"
#include "Collimator.h"
#include "BeamData.h"
#include "Material.h"
#include "MaterialDatabase.h"
#include "LatticeFunctions.h"

using namespace std;
/**
 * Collimator database, used to load and store collimator info
 */
class CollimatorDatabase
{

public:

	/**
	 * Constructor
	 */
	CollimatorDatabase(string, MaterialDatabase*, bool use_sigma);

	/**
	 * Collimator setting structure.
	 */
	struct CollimatorData
	{
		string name;            ///Collimator name
		double x_gap;           ///Collimator x-gap
		double y_gap;           ///Collimator y-gap
		double tilt;            ///Collimator tilt
		double x_offset;        ///Collimator x offset
		double y_offset;        ///Collimator y offset
		double j1_tilt;         ///Collimator jaw 1 tilt
		double j2_tilt;         ///Collimator jaw 2 tilt
		double length;          ///Collimator length (m)
		Material* JawMaterial;  ///Collimator material
		double sigma_x;         ///Jaw x opening in number of sigmas - note that this is not sigma_x it is the sigma in the collimation plane
		double sigma_y;         ///Jaw y opening in number of sigmas - note that this is not sigma_y it is the sigma in the plane orthogonal to the collimation plane
		double beta_x;          ///Calculated x beta function at the collimator entrance
		double beta_y;          ///Calculated y beta function at the collimator entrance
		double position;        ///Length along the lattice, used to calculate the beta functions

	};

	/**
	 * Struct for holding FLUKA data created by A. Valloni + HR
	 */
	struct FlukaData
	{
		int id_coll;            ///Collimator ID
		string name;            ///Collimator name
		double position;        ///Length along the lattice, used to calculate the beta functions
		double angle;           ///Collimator angle = tilt [rad]
		double beta_x;          ///Calculated x beta function at the collimator entrance [m]
		double beta_y;          ///Calculated y beta function at the collimator entrance [m]
		double half_gap;        ///Collimator half gap [m]
		string material;        ///Collimator material symbol
		double length;          ///Collimator length [m]
		double sig_x;           ///Beam sigma_x value [m]
		double sig_y;           ///Beam sigma_y value [m]
		double j1_tilt;         ///Collimator jaw 1 tilt [rad]
		double j2_tilt;         ///Collimator jaw 2 tilt [rad]
		double n_sig;           ///Collimator halfgap in sigma

	};

	CollimatorData* CollData;
	size_t number_collimators;
	bool use_sigma;

	double ConfigureCollimators(AcceleratorModel* model, double emittance_x, double emittance_y,
		LatticeFunctionTable* twiss);
	void ConfigureCollimators(AcceleratorModel* model);
//	void SetOneSideTCDQA(AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* TWISS);

	/**
	 * What impact factor do we wish our halo to impact the collimator jaw (in m).
	 * This just enables a sigma value for bunch generation to be calculated
	 */
	void SelectImpactFactor(string pcoll, double impact);

	/**
	 * Set the stream for the collimator settings log.
	 */

	void SetLogFile(ostream& os);

	/**
	 * Enable/disable logging
	 */
	void EnableLogging(bool);

	/**
	 * Set the stream for the collimator errors log.
	 */
	void SetErrorLogFile(ostream& os);

	/**
	 * Enable/disable logging
	 */
	void EnableErrorLogging(bool);

	/**
	 * Do we match the collimator positions at entry AND exit to the beam envelope? (beta function)
	 */
	void MatchBeamEnvelope(bool);

	/**
	 * Enable/disable the collimators to be matched to the reference orbit and crossing angle
	 */
	void MatchReferenceOrbit(bool);

	/**
	 * Enable jaw flatness errors.
	 */
	void EnableJawFlattnessErrors(bool);

	/**
	 * Enable jaw alignment errors.
	 */
	void EnableJawAlignmentErrors(bool);

	// Set one side jaw for TCDQA collimator
	// void EnableOneSideJawTCDQA(bool);

	/**
	 * Set jaw position offset sigma.
	 */
	void SetJawPositionError(double);

	/**
	 * Set jaw position angle sigma.
	 */
	void SetJawAngleError(double);

	/**
	 * Vector to store FlukaData
	 */
	vector<FlukaData*> StoredFlukaData;

	/**
	 * Function to output FlukaDatabase file
	 */
	void OutputFlukaDatabase(std::ostream* os);

protected:

	string PrimaryCollimator;           /// name of collimator where first impact will occur
	double RequestedImpactFactor;       /// Impact factor in m
	double ImpactSigma;                 /// Impact factor at collimator in number of sigmas
	ostream* log;
	bool logFlag;

	ostream* ErrorLog;
	bool ErrorLogFlag;

	/**
	 * Flag for one side TCDQA jaw
	 */
	bool OneSideJawTCDQA;

	// Match the collimator jaws to the beam envelope? - reference orbit (including crossing angles) and beta functions.
	bool EnableMatchBeamEnvelope;
	bool EnableMatchReferenceOrbit;
	bool JawFlattnessErrors;
	bool JawAlignmentErrors;

	/**
	 * Do we enable resistive collimator wakefields?
	 */
	bool EnableResistiveCollimatorWakes;

	double AngleError;
	double PositionError;
private:

};

#endif
