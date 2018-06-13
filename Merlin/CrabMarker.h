/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CrabMarker_h
#define CrabMarker_h 1

#include "merlin_config.h"

#include "Drift.h"

class ComponentTracker;

/**
 * A CrabMarker is used to initiate the CCFailureProcess
 */

class CrabMarker: public Drift
{
public:

	CrabMarker(const string& id, double len);
	//Overloaded constructor takes the phase advances in x and y
	CrabMarker(const string& id, double len, double mux, double muy);

	/**
	 *	Returns the type string for this component.
	 */
	virtual const string& GetType() const;

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 */
	virtual int GetIndex() const;

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Rotates the component 180 degrees about its local Y axis.
	 */
	virtual void RotateY180();

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	/**
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;

	double GetMuX() const
	{
		return mu_x;
	}
	double GetMuY() const
	{
		return mu_y;
	}
	void SetMuX(double mux)
	{
		mu_x = mux;
	}
	void SetMuY(double muy)
	{
		mu_y = muy;
	}

private:
	double mu_x;
	double mu_y;

};

#endif
