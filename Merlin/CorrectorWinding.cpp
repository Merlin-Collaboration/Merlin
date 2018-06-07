/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "CorrectorWinding.h"

CorrectorWinding::CorrectorWinding(RectMultipole& aMagnet) :
	ModelElement(aMagnet.GetName()), magnet(&aMagnet)
{
}

void CorrectorWinding::SetBx(double value)
{
	if((*magnet).GetField().GetFieldScale() == 0)
	{
		return;
	}

	double l = magnet->GetLength();
	(*magnet).GetField().SetComponent(0, 0, value / l);
}

void CorrectorWinding::SetBy(double value)
{
	if((*magnet).GetField().GetFieldScale() == 0)
	{
		return;
	}

	double l = magnet->GetLength();
	(*magnet).GetField().SetComponent(0, value / l, 0);
}

double CorrectorWinding::GetBx() const
{
	return (*magnet).GetField().GetComponent(0).imag() * magnet->GetLength();
}

double CorrectorWinding::GetBy() const
{
	return (*magnet).GetField().GetComponent(0).real() * magnet->GetLength();
}

const string& CorrectorWinding::GetType() const
{
	_TYPESTR(CorrectorWinding)
}

ModelElement* CorrectorWinding::Copy() const
{
	return new CorrectorWinding(*this);
}
