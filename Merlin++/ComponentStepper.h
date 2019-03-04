/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ComponentStepper_h
#define ComponentStepper_h 1

#include "merlin_config.h"

class AcceleratorComponent;

/**	Utility class used to calculate the required steps
 * through a specific component.
 */
class ComponentStepper
{
public:
	virtual void SetComponent(AcceleratorComponent& cmp) = 0;

	/// Increments the step distance and returns true if on a step boundary, otherwise false.
	virtual bool Increment(double ds) = 0;

	/// Returns the distance to the next step boundary.
	virtual double DistanceToStepBoundary() const = 0;

};

/** A class of stepper which divides a non-zero length
 * component into number of equal steps.
 */
class ComponentDivider: public ComponentStepper
{
public:
	/** Constructor taking the number of steps to take per
	 * component, together with the minimum step distance
	 * (default 0). If for a given component, the calculate
	 * step length is less than min_step, then the step number
	 * is adjusted to give the closest number of integral steps
	 * to match min_step.
	 */
	ComponentDivider(int ns, double min_step = 0);

	virtual void SetComponent(AcceleratorComponent& cmp);

	/// Increments the step distance and returns true if on a step boundary, otherwise false.
	virtual bool Increment(double ds);

	/// Returns the distance to the next step boundary.
	virtual double DistanceToStepBoundary() const;

private:
	// Data Members for Class Attributes

	/// Current component length
	double s;

	/// Integrated length
	double next_s;
	double delta_s;
	double minStep;
	int nstep;
};

#endif
