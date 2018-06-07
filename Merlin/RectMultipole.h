/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RectMultipole_h
#define RectMultipole_h 1

#include "merlin_config.h"
#include "TemplateComponents.h"
#include "MultipoleField.h"
#include "RectMultipoleField.h"
#include "RectangularGeometry.h"

class ComponentTracker;

/**
 *	Abstract base class for all standard rectangular
 *	rectangular multipole magnets that are typically found
 *	in accelerator systems.
 */

class RectMultipole: public RectMultipoleField
{
public:

	/**
	 *	Sets the primary field strength for this multipole.
	 *	Units are Tesla/meter^n, where n is the primary pole
	 *	type.
	 */
	void SetFieldStrength(double b);

	/**
	 *	Returns the primary field strength for this multipole.
	 *	Units are Tesla/meter^n, where n is the primary pole
	 *	type.
	 *	@return Primary field strength
	 */
	double GetFieldStrength() const;

	/**
	 *	Returns the primary pole number for this multipole.
	 *	@return Primary pole number
	 */
	int GetPrimaryPoleNo() const;

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
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

	static const int ID;

protected:

	/**
	 *	Constructor taking the id and the length of the magnet,
	 *	and the definition of a single multipole component. b is
	 *	the field (in Tesla) at the specified radius r0. if
	 *	skew==true, then a skew multipole is constructed
	 *	(default=false).
	 */
	RectMultipole(const string& id, double length, int npole, double b, double r0, bool skew = false);

	/**
	 *	Constructor taking the id and the length of the magnet,
	 *	and the definition of a single multipole component.
	 *	Here, b is the n-th field derivative  (in Tesla/meter^n,
	 *	where n is the pole number), evaluated at r=0.  if
	 *	skew==true, then a skew multipole is constructed
	 *	(default=false).
	 */
	RectMultipole(const string& id, double len, int np, double b, bool skew = false);

private:
	/**
	 *	The primary pole for this multipole.
	 */
	int np;
};

inline void RectMultipole::SetFieldStrength(double b)
{
	GetField().SetFieldScale(b);
}

inline double RectMultipole::GetFieldStrength() const
{
	return GetField().GetFieldScale();
}

inline int RectMultipole::GetPrimaryPoleNo() const
{
	return np;
}

#endif
