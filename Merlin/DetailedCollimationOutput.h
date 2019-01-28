/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef DetailedCollimationOutput_h
#define DetailedCollimationOutput_h 1

#include <string>
#include <vector>

#include "CollimationOutput.h"
#include "AcceleratorComponent.h"
#include "PSTypes.h"
#include "StringPattern.h"

namespace ParticleTracking
{

/**
 * DetailedCollimationOutput is a detailed CollimationOutput
 *
 * DetailedCollimationOutput outputs each individual particle loss, rather
 * than binning by position. The output includes the previous scatter
 * location, so it is useful to diagnose the cause of any given loss.
 *
 */
class DetailedCollimationOutput: public CollimationOutput
{
public:

	DetailedCollimationOutput();
	~DetailedCollimationOutput();

	/**
	 * Not needed for DetailedCollimationOutput, as no binning occurs.
	 */
	virtual void Finalise()
	{
	}
	virtual void Output(std::ostream* os);
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

	/**
	 * Add an element name to record at.
	 *
	 * Can be a pattern, e.g:
	 *
	 *     AddIdentifier("*TCP*")
	 *
	 * to record at elements with "TCP" in their
	 * name.
	 */
	virtual void AddIdentifier(const std::string e);

private:
	std::vector<StringPattern> ids;

};

}

#endif
