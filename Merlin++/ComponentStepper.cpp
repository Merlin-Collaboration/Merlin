/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "utils.h"
#include "AcceleratorComponent.h"
#include "ComponentStepper.h"

ComponentDivider::ComponentDivider(int ns, double min_step) :
	s(0), next_s(0), delta_s(0), minStep(min_step), nstep(ns)
{
}

void ComponentDivider::SetComponent(AcceleratorComponent& cmp)
{
	double l = cmp.GetLength();
	if(l != 0)
	{
		delta_s = l / nstep;
		if(delta_s < minStep)
		{
			int ns = (l / minStep) + 1;
			delta_s = l / ns;
		}
		next_s = delta_s;
	}
	else
	{
		next_s = 0;
	}
	s = 0;
}

bool ComponentDivider::Increment(double ds)
{
	s += ds;
	if(fequal(s, next_s))
	{
		next_s += delta_s;
		return true;
	}
	return false;
}

double ComponentDivider::DistanceToStepBoundary() const
{
	return next_s - s;
}
