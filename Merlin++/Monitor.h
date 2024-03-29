/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Monitor_h
#define Monitor_h 1

#include "merlin_config.h"
#include "SimpleDrift.h"

class Bunch;
class ComponentTracker;

/**
 *	An arbitrary monitor. A monitor can represent any
 *	diagnostic which is typically found in an accelerator.
 *	Although objects of class Monitor can be instantiated to
 *	represent generic monitors, a virtual function Make
 *	Measurement() is provided which can be overridden by
 *	derived Monitor classes which represent more specific
 *	diagnostic types.
 */

class Monitor: public SimpleDrift
{
public:

	Monitor();

	/**
	 *	Constructor taking the identifier for the diagnostic,
	 *	together with its length (default=0) and measurement
	 *	point (default = 0, mid-point).
	 */
	Monitor(const string& id, double len = 0, double mpt = 0);

	virtual ~Monitor();

	/**
	 *	Pure virtual function. Makes a measurement on the
	 *	supplied Beam object. Concrete diagnostics must supply
	 *	this function.
	 */
	virtual void MakeMeasurement(const Bunch&);

	/**
	 *	Sets the position of the measurement point on the local
	 *	geometry. Must be in the range -length/2 to +length/2.
	 */
	void SetMeasurementPt(double mpt);

	/**
	 *	Return the current measurement point.
	 */
	double GetMeasurementPt() const;

	/**
	 *	Returns true if this diagnostic's state is active.
	 */
	bool IsActive() const;

	/**
	 *	Rotate the diagnostic 180 degrees about the local Y axis.
	 */
	virtual void RotateY180();

	/**
	 *	Returns true if the diagnostic is rotated (i.e. the x
	 *	measurement axis is reflected).
	 */
	bool IsReflected() const;

	/**
	 *	Sets the rotated status for this diagnostic to rflg.
	 *	Returns the previous state.
	 */
	bool SetReflected(bool rflg);

	/**
	 *	Sets the local  active state of this monitor. Returns
	 *	the previous state. Note that this method only modifies
	 *	the state of this monitor, and has no affect on the all_
	 *	inactive flag.
	 */
	bool SetActive(bool b);

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string for this component.
	 */
	virtual const string& GetType() const;

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	/**
	 *	Used to turn all monitors off.
	 */
	static bool all_inactive;

	static const int ID;

private:

	/**
	 *	The measurement point of the monitor. Refers to the
	 *	exact point at which the data is to be recorded.
	 */
	double mp;

	/**
	 *	True if the x measurement axis is reflected.
	 *
	 *	@retval true If x measurement axis is reflected
	 *	@retval False If x measurement axis is not reflected
	 */
	bool reflected;

	bool active;
};

inline bool Monitor::IsActive() const
{
	return active && !all_inactive;
}

inline bool Monitor::IsReflected() const
{
	return reflected;
}

inline bool Monitor::SetReflected(bool rflg)
{
	bool tmp = reflected;
	reflected = rflg;
	return tmp;
}

inline bool Monitor::SetActive(bool b)
{
	bool tmp = active;
	active = b;
	return tmp;
}

#endif
