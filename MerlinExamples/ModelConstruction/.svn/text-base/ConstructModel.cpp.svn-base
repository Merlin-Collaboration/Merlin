// ConstructModel.cpp
//
// Merlin 3.0 example code
// v1.0 Nick Walker 21.10.2004
//-------------------------------------------------------------------------


// An example of constructing a Merlin model using the MADInterface class.
//-------------------------------------------------------------------------

// This example uses the MADInterface object to read and construct a Merlin
// model. The source file for the lattice is generated from MAD using the
// MAD OPTICS command. The format of the command is as follows:
//
// SELECT, OPTICS, CLEAR 
// SELECT, OPTICS, #S/#E 
// OPTICS, COLUMNS= NAME, KEYWORD, S, L, K0L, E1, E2, & 
//         K1L, K2L, K3L, TILT, TYPE, & 
//         FILENAME= "output_fname" 
//
// Note that MAD elements that are not supported by Merlin 3.0 are treated
// as drifts by default (a warning message is printed to cout when one is
// encountered).

#include "ConstructModel.h"
#include "MADInterface/MADInterface.h"
#include <iostream>
//#include <fstream>

AcceleratorModel* ConstructModel(const string& fname, double energyInGeV, const string& logFile)
{
	// Construct the parser
	MADInterface madi(fname,energyInGeV);

	// Open a file stream for the log file
	// (useful for debugging)
	std::ofstream logOS(logFile.c_str());
	madi.SetLogFile(logOS);
	madi.SetLoggingOn();

	// The following methods tell MADInterface how to deal
	// with some special cases

	madi.HonourMadStructure(false);
	// This tells MADInterface to ignore the MAD nested BEAMLINE
	// structure, with the exception of BEAMLINES with labels
	// prefixed with X_, where one of the following special
	// SequenceFrame objects are used:
	//    M     MechanicalMover object
	//    S     SimpleMount
	//    G     GirderMount
    // Note the X_ prefix is stripped from the label.

	madi.IgnoreZeroLengthType ("MARKER");
	// In this example we will ignore zero-length MARKER elements

	// madi.TreatTypeAsDrift("RFCAVITY");
	// In this example we treat cavities as simple drifts.

	return madi.ConstructModel();
}



