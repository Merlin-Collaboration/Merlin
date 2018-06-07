/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "SectorBend.h"
#include "ComponentTracker.h"

#include "PhysicalConstants.h"
#include <utility>

/**
 * \class SectorBend::PoleFace
 */
SectorBend::PoleFace::PoleFace(double angle, double f_int, double hg, double face_type) :
	rot(angle), fint(f_int), hgap(hg), type(face_type)
{
}

/**
 * \class SectorBend::PoleFaceInfo
 */
void SectorBend::PoleFaceInfo::Copy(const PoleFaceInfo& rhs)
{
	entrance = rhs.entrance ? new PoleFace(*(rhs.entrance)) : nullptr;
	if(rhs.entrance == rhs.exit)
	{
		exit = entrance;
	}
	else
	{
		exit = rhs.exit ? new PoleFace(*(rhs.exit)) : nullptr;
	}
}

/**
 * \class SectorBend
 */
const int SectorBend::ID = UniqueIndex();

SectorBend::SectorBend(const string& id, double len, double h, double b0) :
	TAccCompGF<ArcGeometry, MultipoleField>(id, new ArcGeometry(len, h), new MultipoleField(0, b0, 1.0, false))
{
}

int SectorBend::GetIndex() const
{
	return ID;
}

const string& SectorBend::GetType() const
{
	_TYPESTR(SectorBend);
}

void SectorBend::RotateY180()
{
	GetField().RotateY180();
	double h = GetGeometry().GetCurvature();
	GetGeometry().SetCurvature(-h);
	std::swap(pfInfo.entrance, pfInfo.exit);
}

void SectorBend::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent);
}

ModelElement* SectorBend::Copy() const
{
	return new SectorBend(*this);
}

void SectorBend::SetPoleFaceInfo(PoleFace* entr, PoleFace* exit)
{
	pfInfo.SetInfo(entr, exit);
}

void SectorBend::SetPoleFaceInfo(PoleFace* pf)
{
	pfInfo.SetInfo(pf);
}

double SectorBend::GetMatchedMomentum(double q) const
{
	using namespace PhysicalConstants;
	double h = GetGeometry().GetCurvature();
	double By = GetField().GetCoefficient(0).real();
	By *= GetField().GetFieldScale();
	return 1.0e-09 * SpeedOfLight * By / h / q;
}

double SectorBend::GetB1() const
{
	return GetField().GetComponent(1).real();
}

void SectorBend::SetB1(double b1)
{
	GetField().SetComponent(1, b1);
}
