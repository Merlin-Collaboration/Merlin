#ifndef _COLLIMATOR_DATABASE_H_
#define _COLLIMATOR_DATABASE_H_

#include <fstream>
#include <string>

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"

#include "BeamModel/BeamData.h"

#include "Collimators/Material.h"
#include "Collimators/MaterialDatabase.h"

#include "RingDynamics/LatticeFunctions.h"

using namespace std;
//Collimator database, used to load and store collimator info
class CollimatorDatabase
{

public:

//Constructor
CollimatorDatabase(string, MaterialDatabase*, bool use_sigma);


//Collimator setting structure.
struct CollimatorData
{
	string name;			//Collimator name
	double x_gap;			//Collimator x-gap
	double y_gap;			//Collimator y-gap
	double tilt;			//Collimator tilt
	double x_offset;		//Collimator x offset
	double y_offset;		//Collimator y offset
	double j1_tilt;			//Collimator jaw 1 tilt
	double j2_tilt;			//Collimator jaw 2 tilt
	double length;			//Collimator length (m)
	Material* JawMaterial;	//Collimator material
	double sigma_x;			//Jaw x opening in number of sigmas
	double sigma_y;			//Jaw y opening in number of sigmas
	double beta_x;			//Calculated x beta function at the collimator entrance
	double beta_y;			//Calculated y beta function at the collimator entrance
	double position;		//Length along the lattice, used to calculate the beta functions
};

	CollimatorData* CollData;
	size_t number_collimators;
	bool use_sigma;



	double ConfigureCollimators(AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss);
	void ConfigureCollimators(AcceleratorModel* model);
//	void SetOneSideTCDQA(AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss);

	//What impact factor do we wish our halo to impact the collimator jaw (in m).
	//This just enables a sigma value for bunch generation to be calculated
	void SelectImpactFactor(string pcoll, double impact);

	//Set the stream for the collimator settings log.
	void SetLogFile (ostream& os);

	//Enable/disable logging
	void EnableLogging(bool);

	//Set the stream for the collimator errors log.
	void SetErrorLogFile (ostream& os);

	//Enable/disable logging
	void EnableErrorLogging(bool);

	//Do we match the collimator positions at entry AND exit to the beam envelope? (beta function)
	void MatchBeamEnvelope(bool);

	//Enable/disable the collimators to be matched to the reference orbit and crossing angle
	void MatchReferenceOrbit(bool);

	//Enable jaw flatness errors.
	void EnableJawFlattnessErrors(bool);

	//Enable jaw alignment errors.
	void EnableJawAlignmentErrors(bool);

	//Set one side jaw for TCDQA collimator
	//void EnableOneSideJawTCDQA(bool);

	//Set jaw position offset sigma.
	void SetJawPositionError(double);

	//Set jaw position angle sigma.
	void SetJawAngleError(double);
	
	
protected:

	string PrimaryCollimator;			// name of collimator where first impact will occur
	double RequestedImpactFactor;		// Impact factor in m
	double ImpactSigma;					// Impact factor at collimator in number of sigmas
	ostream* log;
	bool logFlag;

	ostream* ErrorLog;
	bool ErrorLogFlag;

	// Flag for one side TCDQA jaw
	bool OneSideJawTCDQA; 
	
	// Match the collimator jaws to the beam envelope? - reference orbit (including crossing angles) and beta functions.
	bool EnableMatchBeamEnvelope;
	bool EnableMatchReferenceOrbit;
	bool JawFlattnessErrors;
	bool JawAlignmentErrors;

	// Do we enable resistive collimator wakefields?
	bool EnableResistiveCollimatorWakes;
	
	double AngleError;
	double PositionError;
private:

};

#endif
