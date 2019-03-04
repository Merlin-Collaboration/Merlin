/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ParticleMapComponent_h
#define ParticleMapComponent_h 1

#include "merlin_config.h"
#include "AcceleratorComponent.h"
#include "ParticleMap.h"

namespace ParticleTracking
{

class ParticleBunch;

/**
 *	A special AcceleratorComponent that allows an arbitrary
 *	map (ParticleMap) to be placed into the accelerator
 *	model.
 */
class ParticleMapComponent: public AcceleratorComponent
{
public:

	ParticleMapComponent(const std::string& id, ParticleMap* pmap, double intB2ds = 0);

	/**
	 *	Return the type string for the element.
	 *	@return Element type string
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 *
	 *	@return Unique index of accelerator component class
	 */
	virtual int GetIndex() const;

	/**
	 *	Rotates the component 180 degrees about its local Y axis.
	 */
	virtual void RotateY180();

	ParticleBunch& Apply(ParticleBunch& bunch) const;

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 *	Returns the integral of B^2 for synchrotron radiation
	 *	applications
	 *
	 *	@return The integral \f$ \int B^2 \mathrm{d} s \f$ for synchrotron
	 *	radiation applications
	 */
	double GetIntB2ds() const;

	/**
	 * Data Members for Class Attributes
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;

private:
	/**
	 * Data Members for Associations
	 */
	ParticleMap* itsMap;

	/**
	 * Data Members for Class Attributes
	 */
	double ib2;

}; // Class ParticleMapComponent

inline ParticleBunch& ParticleMapComponent::Apply(ParticleBunch& bunch) const
{
	return itsMap->Apply(bunch);
}

inline double ParticleMapComponent::GetIntB2ds() const
{
	return ib2;
}

} //end namespace ParticleTracking
#endif
