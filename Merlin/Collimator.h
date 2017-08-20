
/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
*
* Class library version 2.0 (2000)
*
* file Merlin\AcceleratorModel\StdComponent\Collimator.h
* last modified 04/04/01 15:25:43
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for
* Charge Particle Accelerator Simulations
*
* Copyright (c) 2000 by The Merlin Collaboration.
* ALL RIGHTS RESERVED.
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
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

class Collimator : public Drift
{
public:

	// Overloaded constructor
	Collimator (const string& id, double len);
	Collimator (const string& id, double len, double radLength);
	Collimator (const string& id, double len, Material* pp, double P0);
	//~ Collimator (const string& id, double len, Material* pp, Collimation::ScatteringModel* s, double P0);

	/**
	* Returns the length of the collimator in units of its
	* radiation length
	* @return Collimator length
	*/
	double GetNumRadLengths() const
	{
		return GetLength()/Xr;
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
	virtual const string& GetType () const;

	/**
	*	Returns the unique index for this class of accelerator
	*	components.
	*	@return Unique index for the class of accelerator components
	*/
	virtual int GetIndex () const;

	/**
	*	Primary tracking interface. Prepares the specified
	*	Tracker object for tracking this component.
	*/
	virtual void PrepareTracker (ComponentTracker& aTracker);

	/**
	*	Rotates the component 180 degrees about its local Y axis.
	*/
	virtual void RotateY180 ();

	/**
	*	Virtual constructor.
	*/
	virtual ModelElement* Copy () const;

	/**
	*	Unique index for an Accelerator component.
	*/
	static const int ID;

	/**
	* Set collimator ID for Fluka output
	*/
	virtual void SetCollID (int n)
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

	bool scatter_at_this_collimator;

	/**
	* Collimator material
	*/
	Material* p;
	virtual void SetMaterial(Material* pp)
	{
		p = pp;
	};

	// ScatteringModel contains the relevent ScatteringProcess to use when performing scattering
	//~ Collimation::ScatteringModel* scatter;
	//~ virtual void SetScatteringModel(Collimation::ScatteringModel* s);


private:
	double Xr;
};

#endif
