/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CorrectorDipoles_h
#define CorrectorDipoles_h 1

#include "merlin_config.h"
#include "RectMultipole.h"

class ComponentTracker;

/**
 *	A horizontal corrector dipole.
 */

class XCor: public RectMultipole
{
public:

	/**
	 *	Constructor taking the identifier for the corrector, the
	 *	length and the field (in Tesla).
	 */
	XCor(const string& id, double len, double B = 0);

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

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

	static const int ID;
};

/**
 *	A vertical corrector dipole.
 */

class YCor: public RectMultipole
{
public:

	/**
	 *	Constructor taking the identifier for the corrector, the
	 *	length and the field (in Tesla).
	 */
	YCor(const string& id, double len, double B = 0);

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

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

	static const int ID;
};

inline XCor::XCor(const string& id, double len, double B) :
	RectMultipole(id, len, 0, B)
{
}

inline YCor::YCor(const string& id, double len, double B) :
	RectMultipole(id, len, 0, B, true)
{
}

#endif
