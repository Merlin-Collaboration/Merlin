/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/StdComponent/RFStructure.h"

#define _FIELD static_cast<RFAcceleratingField*>(this->itsField)

double RFStructure::GetFrequency () const
{
    return _FIELD->GetFrequency();
}

double RFStructure::GetPhase () const
{
    return _FIELD->GetPhase();
}

void RFStructure::SetAmplitude (double Epk)
{
    _FIELD->SetAmplitude(Epk);
}

double RFStructure::GetAmplitude () const
{
    return _FIELD->GetAmplitude();
}

double RFStructure::GetWavelength () const
{
    return _FIELD->GetWavelength();
}

double RFStructure::GetK () const
{
    return _FIELD->GetK();
}

void RFStructure::SetFrequency (double f)
{
    _FIELD->SetFrequency(f);
}

void RFStructure::SetPhase (double phase)
{
    _FIELD->SetPhase(phase);
}

void RFStructure::SetWavelength (double lambda)
{
    using PhysicalConstants::SpeedOfLight;
    _FIELD->SetFrequency(SpeedOfLight/lambda);
}

void RFStructure::SetK (double k)
{
    using PhysicalConstants::SpeedOfLight;
    _FIELD->SetFrequency(SpeedOfLight*k/twoPi);
}

RFStructure::RFStructure (const string& id, double len, RFAcceleratingField* aField)
        : TAccCompGF_NC<RectangularGeometry,RFAcceleratingField>(id,new RectangularGeometry(len),aField)
{}

double RFStructure::GetBeamVoltage () const
{
    return GetAmplitude()*cos(GetPhase())*GetLength();
}

double RFStructure::GetVoltage() const
{
	return GetAmplitude()*GetLength();
}

void RFStructure::SetVoltage(double v)
{
	SetAmplitude(v/GetLength());
}

Complex RFStructure::GetVoltagePhasor() const
{
	double phi=GetPhase();
	return GetVoltage()*Complex(cos(phi),sin(phi));
}

void RFStructure::SetVoltagePhasor(const Complex& z)
{
	SetVoltage(abs(z));
	SetPhase(arg(z));
}
