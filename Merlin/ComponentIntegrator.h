/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ComponentIntegrator_h
#define ComponentIntegrator_h 1

#include "merlin_config.h"
#include "utils.h"
#include "Range.h"
#include "AcceleratorComponent.h"
#include "utils.h"

/**
 *  Class ComponentIntegrator
 *
 *	Abstract class used by ComponentTracker to track through an
 *	AcceleratorComponent. Integrators encapsulate the
 *	various algorithms used for tracking.
 */

class ComponentIntegrator
{
public:

	ComponentIntegrator();
	virtual ~ComponentIntegrator();

	/**
	 *	Track the entire component.
	 */
	void TrackAll();

	/**
	 *	Tracks a single step ds through the current component.
	 *  Returns the remaining length.
	 *  @return Remaining length of current component
	 */
	virtual double Track(double ds);

	/**
	 *	Returns the component index for this integrator.
	 *	@return Component index for this integrator
	 */
	virtual int GetComponentIndex() const = 0;

	/**
	 *	Returns the total integrated length through the current
	 *	component.
	 *	@return Total integrated length through current component
	 */
	double GetIntegratedLength() const;

	/**
	 *	Returns the remaining length to integrate.
	 *	@return Length remaining to integrate
	 */
	double GetRemainingLength() const;

	/**
	 *	Returns true if the step ds is valid (allowed)
	 *	@retval true If the step ds is valid
	 *	@retval false If the step ds is invalid
	 */
	bool IsValidStep(double ds) const;

	/**
	 *  Returns the remaining valid step range for this integrator.
	 *  Default behaviour returns [0,GetRemainingLength()]
	 *  @return Remaining valid step range for the integrator
	 */
	virtual FloatRange GetValidStepRange() const;

	/**
	 *	Component boundary checks
	 *  Returns true if GetIntegratedLength()==0
	 *  @retval true If `GetIntegratedLength()==0`
	 */
	bool AtEntrance() const;

	/**
	 *  Returns true if step == GetRemainingLength()
	 *  @retval true If `step == GetRemainingLength()`
	 */
	bool AtExit(double step = 0) const;

protected:

	friend class ComponentTracker;

	/**
	 *  Sets the current component to track.
	 *  Can only be called by ComponentTracker via friendship
	 */
	virtual void SetCurrentComponent(AcceleratorComponent& aComponent);

	/**
	 *  Perform the tracking for step ds. Concrete integrators
	 *  must supply this function.
	 */
	virtual void TrackStep(double ds) = 0;

	/**
	 * Function applies entrance field map
	 */
	virtual void TrackEntrance()
	{
	}

	/**
	 * Function applies exit field map
	 */

	virtual void TrackExit()
	{
	}

private:

	double cur_S; /// current integrated length
	double tot_S; /// total length

	/**
	 *  The current component
	 */
	AcceleratorComponent* component;

	// Copy protection
	ComponentIntegrator(const ComponentIntegrator& rhs);
	ComponentIntegrator& operator=(const ComponentIntegrator& rhs);
};

inline ComponentIntegrator::ComponentIntegrator() :
	cur_S(0), tot_S(0), component(nullptr)
{
}

inline ComponentIntegrator::~ComponentIntegrator()
{
}

inline void ComponentIntegrator::TrackAll()
{
	Track(component->GetLength());
}

inline void ComponentIntegrator::SetCurrentComponent(AcceleratorComponent& aComponent)
{
	component = &aComponent;
	cur_S = 0;
	tot_S = component->GetLength();
}

inline double ComponentIntegrator::GetIntegratedLength() const
{
	return cur_S;
}

inline double ComponentIntegrator::GetRemainingLength() const
{
	double s = tot_S - cur_S;
	return fequal(s, 0) ? 0 : s;
}

inline bool ComponentIntegrator::IsValidStep(double ds) const
{
	return GetValidStepRange()(ds);
}

inline double ComponentIntegrator::Track(double ds)
{
	if(AtEntrance())
	{
		TrackEntrance();
	}

	TrackStep(ds);
	cur_S += ds;

	if(AtExit())
	{
		TrackExit();
	}

	return GetRemainingLength();
}

inline bool ComponentIntegrator::AtEntrance() const
{
	return fequal(cur_S, 0.0);
}

inline bool ComponentIntegrator::AtExit(double s) const
{
	return fequal(s + cur_S, tot_S);
}

inline FloatRange ComponentIntegrator::GetValidStepRange() const
{
	return FloatRange(0, GetRemainingLength());
}

#endif
