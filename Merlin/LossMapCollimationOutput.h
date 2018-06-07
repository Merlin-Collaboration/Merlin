/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef LossMapCollimationOutput_h
#define LossMapCollimationOutput_h 1

#include <string>
#include <vector>

#include "CollimationOutput.h"
#include "AcceleratorComponent.h"
#include "ParticleBunch.h"
#include "PSTypes.h"
namespace ParticleTracking
{

class LossMapCollimationOutput: public CollimationOutput
{

public:

	LossMapCollimationOutput(OutputType otype = tencm);
	~LossMapCollimationOutput();

	/**
	 * Finalise will call any sorting algorithms and perform formatting for final output
	 */
	virtual void Finalise();

	/**
	 * Outputs the loss map data to a specified output stream
	 * @param[out] os The stream to output to.
	 */
	virtual void Output(std::ostream* os);

	/**
	 * Called from CollimateProtonProcess::DoScatter to add a particle to the CollimationOutput
	 */
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

	/**
	 * Sets a warm area of the machine.
	 * @param[in] wr A std::pair that contains the start and end location of a warm region. First contains the start location, and second the end.
	 */
	void SetWarmRegion(std::pair<double, double> wr);

	/**
	 * Clears out any previously defined warm regions.
	 */
	void ClearWarmRegions();

	/**
	 * Gets the vector containing the currently set warm regions of the machine.
	 * @return The std::vector of std::pair<double,double> with each entry containing the start and end locations of warm regions.
	 */
	std::vector<std::pair<double, double> > GetWarmRegions() const;

protected:

	//A vector of std::pair containing the start and end of warm regions of the machine. Can be empty. First contains the start location, and second the end.
	std::vector<std::pair<double, double> > WarmRegions;

private:

};

}
#endif
