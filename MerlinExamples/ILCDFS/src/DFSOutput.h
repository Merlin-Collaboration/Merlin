/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_DFSOutput
#define _h_DFSOutput 1

#include "TrackingSimulation.h"
#include <string>
#include <fstream>

// Used to define output for final (emittance) tracking.
class DFSOutput: public SimulationOutput
{
public:

	DFSOutput() :
		fos(nullptr)
	{
	}

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
