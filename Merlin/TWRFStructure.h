/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef TWRFStructure_h
#define TWRFStructure_h 1

#include "merlin_config.h"
#include "RFStructure.h"

class ComponentTracker;

/**
 *	A travelling wave accelerating structure.
 */
class TWRFStructure: public RFStructure
{
public:

	/**
	 *	Constructor taking the label for this cavity (id), the
	 *	cavity length (len) in meters, the frequency (f) in MHz,
	 *	and gradient (Epk) in MV/m and the phase (phi) in
	 *	radians.
	 */
	TWRFStructure(const string& id, double len, double f, double Epk, double phi = 0);
	TWRFStructure(const TWRFStructure&);

	/**
	 *	Returns the type string "TWRFStructure".
	 *	@return Type string TWRFStructure
	 */
	virtual const string& GetType() const;

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 *
	 *	@return Accelerator component class unique index
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

	static const int ID;
};

#endif
