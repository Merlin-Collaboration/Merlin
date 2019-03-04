/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef TransverseRFStructure_h
#define TransverseRFStructure_h 1

#include "merlin_config.h"
#include "RFStructure.h"
#include "TransRFfield.h"

class ComponentTracker;

/**
 *	A travelling wave transverse deflecting structure.
 */
class TransverseRFStructure: public RFStructure
{
public:

	/**
	 *	Constructor taking the label for this cavity (id), the
	 *	cavity length (len) in meters, the frequency (f) in MHz,
	 *	and gradient (Epk) in MV/m and the phase (phi) in
	 *	radians.
	 */
	TransverseRFStructure(const string& id, double len, double f, double Epk, double phi = 0, double theta = 0);

	/**
	 *  Get Field Orientation
	 */
	double GetFieldOrientation() const;

	/**
	 *  Set Field Orientation
	 */
	void SetFieldOrientation(double);

	/**
	 *	Returns the type string "TransverseRFStructure".
	 *	@return Type string TransverseRFStructrure
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

protected:
private:
	// Private copy constructor
	TransverseRFStructure(const TransverseRFStructure&);
};

inline double TransverseRFStructure::GetFieldOrientation() const
{
	return static_cast<TransverseRFfield*>(itsField)->GetFieldOrientation();
}

inline void TransverseRFStructure::SetFieldOrientation(double t)
{
	static_cast<TransverseRFfield*>(itsField)->SetFieldOrientation(t);
}

#endif
