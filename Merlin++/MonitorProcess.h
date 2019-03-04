/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MonitorProcess_h
#define MonitorProcess_h 1

#include "ParticleBunchProcess.h"
#include <vector>
#include <string>
#include "ParticleBunch.h"

namespace ParticleTracking
{
/**
 * A process for recording particle coordinates at elements
 *
 * Can be attached to the tracker to record particle coordinates at all
 * or specific elements to files.
 */
class MonitorProcess: public ParticleBunchProcess
{
private:
	vector<string> dump_at_elements;
	string file_prefix;
	unsigned int count;

public:
	/**
	 * Create a MonitorProcess
	 *
	 * \param aID Process ID
	 * \param prio Process priority
	 * \param prefix Output file name prefix
	 *
	 */
	MonitorProcess(const string& aID = "MONITOR", int prio = 0, const string& prefix = "");

	/// Set the output file name prefix
	void SetPrefix(const string& prefix);

	/// Add element at which to record
	void AddElement(const string e);
	void InitialiseProcess(Bunch&  bunch);
	void DoProcess(const double ds);
	double GetMaxAllowedStepSize() const;
	void SetCurrentComponent(AcceleratorComponent& component);

};

} // end namespace ParticleTracking
#endif
