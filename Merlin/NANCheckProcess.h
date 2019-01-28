/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef NANCheckProcess_h
#define NANCheckProcess_h 1

#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"
#include <set>

namespace ParticleTracking
{

/**
 * A diagnostic process for catching particles with invalid coordinates
 *
 * Particle can get invalid coordinates due to extreme values of high
 * amplitude particles or bugs in MERLIN or user code. These may eventually
 * cause a crash. This process can be used to catch the invalid (Not a
 * Number) values immediately, so that the cause can be identified.
 *
 * There are several option flags. (disabled by default)
 *
 * detailed: output the previous good coordinates for the particle (requires
 * additional memory to store).
 *
 * cull: removes the invalid particles from the bunch.
 *
 * halt: stops the simulation when an invalid particle is found.
 *
 */
class NANCheckProcess: public ParticleBunchProcess
{
public:

	NANCheckProcess(const string& aID = "NAN Check", int prio = -1);
	void InitialiseProcess(Bunch&  bunch);
	void DoProcess(const double ds);
	double GetMaxAllowedStepSize() const;
	void SetCurrentComponent(AcceleratorComponent& component);

	/// Enable or disable detailed mode
	void SetDetailed(bool enable = true)
	{
		detailed = enable;
	}
	/// Enable or disable cull mode
	void SetCullNAN(bool enable = true)
	{
		cull = enable;
	}
	/// Enable or disable halt mode
	void SetHaltNAN(bool enable = true)
	{
		halt = enable;
	}
private:
	bool detailed;
	bool cull;
	bool halt;

	void Report(int id) const;
	void DoCull();
	PSvectorArray start_coords; // particles at start of element
	PSvectorArray prev_coords; // particles at start of element
	std::set<double> reported;
};

} // end namespace ParticleTracking
#endif
