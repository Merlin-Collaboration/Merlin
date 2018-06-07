/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Marker_h
#define Marker_h 1

#include "merlin_config.h"
#include "AcceleratorComponent.h"

class ComponentTracker;

/**
 *	Special Marker component. This does not represent a
 *	physical component, and is provided for uses wishing to
 *	"mark" specific locations in the lattice. It has no
 *	field or geometry associated with it.
 */

class Marker: public AcceleratorComponent
{
public:

	/**
	 *	Constructor
	 */
	explicit Marker(const std::string& id);

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 *	@return Unique index for the class of accelerator components
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string for this component.
	 *	@return Type string for component
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor. - FIXME
	 */
	virtual ModelElement* Copy() const;

	/**
	 *	Rotates the component 180 degrees about its local Y axis.
	 */
	virtual void RotateY180();

	/**
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;
};

inline Marker::Marker(const std::string& id) :
	AcceleratorComponent(id, nullptr, nullptr)
{
}

#endif
