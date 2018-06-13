/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Drift_h
#define Drift_h 1

#include "merlin_config.h"
#include "TemplateComponents.h"
#include "SimpleDrift.h"
#include "RectangularGeometry.h"

class ComponentTracker;

/**
 *	A simple drift section.
 */

class Drift: public SimpleDrift
{
public:

	explicit Drift(const string& id, double len = 0);

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

	// Data Members for Class Attributes

	/**
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;
};

#endif
