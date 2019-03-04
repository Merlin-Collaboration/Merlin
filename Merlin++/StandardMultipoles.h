/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef StandardMultipoles_h
#define StandardMultipoles_h 1

#include "merlin_config.h"

#include "RectMultipole.h"

class ComponentTracker;

// type strings for template parameter instantiation

/**
 *	A standard quadrupole magnet.
 */
class Quadrupole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	the quadrupole gradient. Units are Tesla/meter.
	 */
	Quadrupole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a specific radius r0 (in meter).
	 */
	Quadrupole(const string& id, double len, double B, double r0);

	/**
	 *	Initialise aTracker with *this.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the unique index for a Quadrupole
	 *	@return Quadrupole unique index
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string "Quadrupole"
	 *	@return Type string "Quadrupole"
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	/**
	 * Data Members for Class Attributes
	 */
	static const int ID;

protected:
private:
private:
};

//	A standard sextupole magnet.

class Sextupole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	a sextupole gradient. Units are Tesla/meter^2.
	 */
	Sextupole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a specific radius r0 (in meter).
	 */
	Sextupole(const string& id, double len, double B, double r0);

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the unique index for a Sextupole
	 *	@return Sextupole unique index
	 */

	virtual int GetIndex() const;

	/**
	 *	Returns the type string "Sextupole".
	 *	@return Type string sextupole
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	// Data Members for Class Attributes

	static const int ID;

protected:
private:
private:
};

/**
 *	 A standard skew-quadrupole magnet.
 */
class SkewQuadrupole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	a quadrupole gradient. Units are Tesla/meter.
	 */
	SkewQuadrupole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a speciific radius r0 (in meter).
	 */
	SkewQuadrupole(const string& id, double len, double B, double r0);

	/**
	 *	Initialise aTracker to track *this.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Return the index for a SkewQuadrupole.
	 *	@return Index for a SkewQuadrupole
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string "SkewQuadrupole".
	 *	@return Type string SkewQuadrupole
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	// Data Members for Class Attributes

	static const int ID;

protected:
private:
private:
};

/**
 *	A standard Octupole magnet.
 */
class Octupole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	a octupole gradient. Units are Tesla/meter^3.
	 */
	Octupole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a specific radius r0 (in meter).
	 */
	Octupole(const string& id, double len, double B, double r0);

	/**
	 *	Initialise aTracker to track *this.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the index for an Octupole.
	 *	@return Index for an Octupole
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string "Octupole".
	 *   @return Type string Octupole
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	// Data Members for Class Attributes

	static const int ID;

protected:
private:
private:
};

class SkewSextupole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	a sextupole gradient. Units are Tesla/meter^2.
	 */
	SkewSextupole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a specific radius r0 (in meter).
	 */
	SkewSextupole(const string& id, double len, double B, double r0);

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the index for a "SkewSextupole".
	 *	@return Index for an SkewSextupole
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string "Skewsextupole".
	 *   @return Type string Skewsextupole
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	// Data Members for Class Attributes

	static const int ID;

protected:
private:
private:
};

/**
 *	A standard (normal) decapole magnet.
 */
class Decapole: public RectMultipole
{
public:
	/**
	 *	Constructor. The primary field component is specified as
	 *	a octupole gradient. Units are Tesla/meter^4.
	 */
	Decapole(const string& id, double len, double dnB);

	/**
	 *	Constructor. The primary field component is specified as
	 *	a field B (in Tesla) at a specific radius r0 (in meter).
	 */
	Decapole(const string& id, double len, double B, double r0);

	/**
	 *	Initialise aTracker to track *this.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the index for a decapole.
	 *	@return Index for a decapole
	 */
	virtual int GetIndex() const;

	/**
	 *	Returns the type string "Decapole".
	 *	@return Type string Decapole
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	static const int ID;

protected:
private:
private:
};

// Class Quadrupole

// Class Sextupole

// Class SkewQuadrupole

// Class Octupole

// Class SkewSextupole

// Class Decapole

#endif
