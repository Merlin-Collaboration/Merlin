/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Collimator_h
#define Collimator_h 1

#include "merlin_config.h"
#include "Drift.h"
#include "Material.h"
#include "ScatteringModel.h"

class ComponentTracker;

/**
 * A collimator represents a scattering element in the beamline. Collimator objects
 * are optically drifts, but when associated with an Aperture, they have the
 * material property of radiation length, which is used by specific scattering
 * process to approximate the interaction of particles with the collimator material.
 */

class Collimator: public Drift
{
public:

	Collimator(const string& id, double len);
	Collimator(const string& id, double len, double radLength);
	Collimator(const string& id, double len, Material* pp, double P0);

	/**
	 * Returns the length of the collimator in units of its radiation length
	 * @return Collimator length
	 */
	double GetNumRadLengths() const
	{
		return GetLength() / Xr;
	}

	/**
	 * Returns the material radiation length (meter)
	 * @return Radiation material length (m)
	 */
	double GetMaterialRadiationLength() const
	{
		return Xr;
	}

	/**
	 *	Returns the type string for this component.
	 *	@return String type of component
	 */
	virtual const string& GetType() const;

	/**
	 *	Returns the unique index for this class of accelerator components.
	 *	@return Unique index for the class of accelerator components
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

	/**
	 * Set collimator ID for Fluka output
	 */
	virtual void SetCollID(int n)
	{
		Coll_ID = n;
	}

	/**
	 * Get collimator ID for Fluka output
	 * @return Collimator ID for Fluka output
	 */
	virtual int GetCollID()
	{
		return Coll_ID;
	}

	/**
	 * Collimator material
	 */
	Material* material;

	virtual void SetMaterial(Material* mat)
	{
		material = mat;
	}

	Material* GetMaterial() const
	{
		return material;
	}

private:
	double Xr;
	bool scatter_at_this_collimator;
};

#endif
