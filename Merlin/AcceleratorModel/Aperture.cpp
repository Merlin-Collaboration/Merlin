#include "AcceleratorModel/Aperture.h"

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
