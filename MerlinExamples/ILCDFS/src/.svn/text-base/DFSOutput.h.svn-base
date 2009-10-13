/////////////////////////////////////////////////////////////////////////
// Class DFSOutput
// Used to define output for final (emittance) tracking.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 18:32:23 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_DFSOutput
#define _h_DFSOutput 1

#include "BeamDynamics/TrackingSimulation.h"
#include <string>
#include <fstream>

class DFSOutput : public SimulationOutput {
public:

	DFSOutput() : fos(0) {}

	// Close any existing file stream and open
	// a new output file.
	void NewFile(const std::string&);

protected:

	virtual void Record(const ComponentFrame* frame, const Bunch* bunch);
    virtual void RecordInitialBunch(const Bunch* bunch);
	virtual void RecordFinalBunch(const Bunch* bunch);

private:
	void Record(const std::string& id, const Bunch* bunch);
	std::ofstream* fos;
};

#endif






