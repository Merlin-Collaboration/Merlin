/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ConstructSrot.h"
#include "PatchFrame.h"

ComponentFrame* ConstructSrot(double angle, const std::string& name)
{
	GeometryPatch* gp = new GeometryPatch;
	gp->RotateZ(angle);
	return new PatchFrame(gp, name);
}

ComponentFrame* ConstructXrot(double angle, const std::string& name)
{
	if(angle)
	{
		GeometryPatch* gp = new GeometryPatch;
		gp->RotateX(angle);
		return new PatchFrame(gp, name);
	}
	else
	{
		return new PatchFrame(nullptr, name);
	}
}
