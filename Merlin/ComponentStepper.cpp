/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 *
 * Class library version 2.0 (1999)
 *
 * file Merlin\BeamDynamics\Utilities\ComponentStepper.cpp
 * last modified 01/07/02 16:36:22
 *
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for
 * Charge Particle Accelerator Simulations
 * Copyright (c) 2001 by The Merlin Collaboration.
 * - ALL RIGHTS RESERVED -
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

#include "utils.h"
#include "AcceleratorComponent.h"
#include "ComponentStepper.h"

ComponentDivider::ComponentDivider (int ns, double min_step)
	:s(0),next_s(0),delta_s(0),minStep(min_step),nstep(ns)
{
}


void ComponentDivider::SetComponent (AcceleratorComponent& cmp)
{
	double l = cmp.GetLength();
	if(l!=0)
	{
		delta_s = l/nstep;
		if(delta_s<minStep)
		{
			int ns = (l/minStep)+1;
			delta_s = l/ns;
		}
		next_s = delta_s;
	}
	else
	{
		next_s = 0;
	}
	s=0;
}

bool ComponentDivider::Increment (double ds)
{
	s+=ds;
	if(fequal(s,next_s))
	{
		next_s+=delta_s;
		return true;
	}
	return false;
}

double ComponentDivider::DistanceToStepBoundary () const
{
	return next_s-s;
}
