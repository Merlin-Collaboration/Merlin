/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Aperture.h"

void Aperture::printout(std::ostream& out) const
{
	out << GetApertureType();
}

Material* Aperture::GetMaterial() const
{
	return ApertureMaterial;
}

void Aperture::SetMaterial(Material* m)
{
	ApertureMaterial = m;
}

inline bool Aperture::PointInside (const Point3D& p) const
{
	return PointInside(p.x,p.y,p.z);
}

std::ostream& operator<< (std::ostream& out, const Aperture& ap)
{
	ap.printout(out);
	return out;
}

Aperture::~Aperture () {}
