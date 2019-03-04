/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ProcessStepManager_h
#define ProcessStepManager_h 1

#include "merlin_config.h"
#include <list>
#include <ostream>

class AcceleratorComponent;
class BunchProcess;
class Bunch;

/**
 * Responsible for coordinating the tracking of a bunch
 * through a single AcceleratorComponent. Tracking occurs
 * by the application of a series of processes, which
 * effectively integrate the bunch motion along the
 * beamline.
 */
class ProcessStepManager
{
public:

	//	Construction/destruction.
	ProcessStepManager();
	~ProcessStepManager();

	/**
	 * Initialise the step manager with the specified (initial) bunch.
	 */
	void Initialise(Bunch& bunch);

	/**
	 * Track the specified component. The current bunch object
	 * is updated accordingly.
	 */
	void Track(AcceleratorComponent& component);

	/**
	 * Returns the total length integrated since the last call to Initialise(Bunch&).
	 * @return Total length integrate since last call to `Initialise(Bunch&)`
	 */
	double GetIntegratedLength();

	/**
	 * Add a process.
	 */
	void AddProcess(BunchProcess* aProcess);

	/**
	 * Remove aProcess from the current process table. Returns
	 * true if aProcess was present, otherwise false.
	 *
	 * @retval true If aProcess present
	 * @retval false If aProcess not present
	 */
	bool RemoveProcess(BunchProcess* aProcess);

	/**
	 * Remove and destroys all processes in the current process table.
	 */
	void ClearProcesses();

	void SetLogStream(std::ostream* os);

private:

	/**
	 * The current integrated length.
	 */
	double total_s;

	std::ostream* log;

	/**
	 * list of processes in order of priority.
	 */
	std::list<BunchProcess*> processTable;

	//Copy protection
	ProcessStepManager(const ProcessStepManager& rhs);
	ProcessStepManager& operator=(const ProcessStepManager& rhs);
};

#endif
