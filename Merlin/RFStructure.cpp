/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "RFStructure.h"

#define _FIELD static_cast<RFAcceleratingField*>(this->itsField)

double RFStructure::GetFrequency() const
{
	return _FIELD->GetFrequency();
}

double RFStructure::GetPhase() const
{
	return _FIELD->GetPhase();
}

void RFStructure::SetAmplitude(double Epk)
{
	_FIELD->SetAmplitude(Epk);
}

double RFStructure::GetAmplitude() const
{
	return _FIELD->GetAmplitude();
}

double RFStructure::GetWavelength() const
{
	return _FIELD->GetWavelength();
}

double RFStructure::GetK() const
{
	return _FIELD->GetK();
}

void RFStructure::SetFrequency(double f)
{
	_FIELD->SetFrequency(f);
}

void RFStructure::SetPhase(double phase)
{
	_FIELD->SetPhase(phase);
}

void RFStructure::SetWavelength(double lambda)
{
	using PhysicalConstants::SpeedOfLight;
	_FIELD->SetFrequency(SpeedOfLight / lambda);
}

void RFStructure::SetK(double k)
{
	using PhysicalConstants::SpeedOfLight;
	_FIELD->SetFrequency(SpeedOfLight * k / twoPi);
}

RFStructure::RFStructure(const string& id, double len, RFAcceleratingField* aField) :
	TAccCompGF_NC<RectangularGeometry, RFAcceleratingField>(id, new RectangularGeometry(len), aField)
{
}

double RFStructure::GetBeamVoltage() const
{
	return GetAmplitude() * cos(GetPhase()) * GetLength();
}

double RFStructure::GetVoltage() const
{
	return GetAmplitude() * GetLength();
}

void RFStructure::SetVoltage(double v)
{
	SetAmplitude(v / GetLength());
}

Complex RFStructure::GetVoltagePhasor() const
{
	double phi = GetPhase();
	return GetVoltage() * Complex(cos(phi), sin(phi));
}

void RFStructure::SetVoltagePhasor(const Complex& z)
{
	SetVoltage(abs(z));
	SetPhase(arg(z));
}
