/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BunchProcess_h
#define BunchProcess_h 1

#include "merlin_config.h"
#include <string>

class AcceleratorComponent;
class Bunch;

using namespace std;

/**
 *	A single process which is applied to a Bunch object
 *	during tracking. A Process can represent any physical
 *	(e.g. particle transport) or abstract mechanism (e.g.
 *	output) which may be applied to a bunch at a specific
 *	AcceleratorComponent. Concrete processes may be generic
 *	to all Bunch model types , or may be specific to a
 *	concrete bunch model (see for example ParticleBunch and
 *	ParticleBunchTransortProcess).
 */

class BunchProcess
{
public:

	/**
	 *	Constructing taking the process ID and its priority
	 *	(default =0, high).
	 */
	explicit BunchProcess(const string& anID, int aPriority = 0);

	virtual ~BunchProcess();

	/**
	 *	Set the priority. Highest priority is 0, followed by
	 *	1,2,...etc.
	 */
	void SetPriority(int p);

	/**
	 *	Get the priority. Highest priority is 0, followed by
	 *	1,2,...etc.
	 */
	int GetPriority() const;

	/**
	 *	Initialise this process with the specified Bunch.
	 */
	virtual void InitialiseProcess(Bunch& bunch) = 0;

	/**
	 *	Sets the current accelerator component. This function
	 *	should be called just before tracking of the component
	 *	begins. Concrete processes should override this function
	 *	to perform component and process dependent
	 *	initialisation.
	 */
	virtual void SetCurrentComponent(AcceleratorComponent& component)
	{
		currentComponent = &component;
	}

	/**
	 *	Preform the process for the specified step ds.
	 */
	virtual void DoProcess(double ds) = 0;

	/**
	 *	Returns the current maximum step length for this process.
	 *	@return Current maximum step length for the process
	 */
	virtual double GetMaxAllowedStepSize() const = 0;

	/**
	 *	Returns true if this process is active.
	 *	@retval true If process is active
	 *	@retval false If process is inactive
	 */
	bool IsActive() const;

	const string& GetID() const;

protected:

	bool active;
	AcceleratorComponent* currentComponent;

private:

	/*
	 * The following two lines disable the copying of BunchProcess objects via the copy constructor or via assignment.
	 */
	BunchProcess(const BunchProcess& bp);
	BunchProcess& operator=(const BunchProcess& bp);

	string ID;
	int priority;
};

/**
 * Template class for defining bunch specific BunchProcess classes
 */

template<class B> class TBunchProc: public BunchProcess
{
public:

	explicit TBunchProc(const std::string& anID, int aPriority = 0) :
		BunchProcess(anID, aPriority)
	{
	}

	/**
	 * Sets the current bunch, if bunch is of type bunch
	 */
	virtual void InitialiseProcess(Bunch& bunch)
	{
		currentBunch = dynamic_cast<B*>(&bunch);
	}

protected:

	B* currentBunch;
};

inline BunchProcess::BunchProcess(const string& anID, int aPriority) :
	active(false), currentComponent(nullptr), ID(anID), priority(aPriority)
{
}

inline BunchProcess::~BunchProcess()
{
}

inline void BunchProcess::SetPriority(int p)
{
	priority = p;
}

inline int BunchProcess::GetPriority() const
{
	return priority;
}

inline bool BunchProcess::IsActive() const
{
	return active;
}

inline const string& BunchProcess::GetID() const
{
	return ID;
}

#endif
