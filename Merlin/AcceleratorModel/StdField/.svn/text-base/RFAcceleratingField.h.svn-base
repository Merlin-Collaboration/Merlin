/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:36:44 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef RFAcceleratingField_h
#define RFAcceleratingField_h 1

#include "merlin_config.h"
#include "NumericalUtils/NumericalConstants.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/EMField.h"

//	An RF accelerating field. An RF accelerating field is
//	characterised by a sinusoidal electric field,
//	which is completely defined by its frequency, amplitude
//	and phase. RFAcceleratingField implements the basic
//	defining attributes for sinusoidal RF fields.

class RFAcceleratingField : public EMField
{
public:

    //	Constructor taking the amplitude (Epk), the frequency
    //	(f) and the phase of the field.
    RFAcceleratingField (double Epk, double f, double phase = 0);

    //	Set the frequency.
    void SetFrequency (double f);

    //	Returns the frequency.
    double GetFrequency () const;

    //	Set the phase.
    void SetPhase (double phase);

    //	Returns the phase.
    double GetPhase () const;

    //	Set the amplitude.
    void SetAmplitude (double Epk);

    //	Returns the amplitude.
    double GetAmplitude () const;

    //	Calculate the Ez component.
    virtual double Ez (double z, double t) const = 0;

    //	Returns the wavelength of the RF (in meter).
    double GetWavelength () const;

    //	Returns the k value (=2pi/wavelength) for the field.
    double GetK () const;

    // Set acceleration mode
    void SetFullAcceleration (bool full_acceln);

    // Get acceleration mode
    bool FullAcceleration () const;

protected:

    //	Frequency (2*pi*f) in radian/s.
    double w;

    //	Peak electric field.
    double E0;

    //	Arbitrary phase angle (radian)
    double phi;

    bool full_acceleration;
};

// Class RFAcceleratingField

inline RFAcceleratingField::RFAcceleratingField (double Epk, double f, double phase)
        : w(twoPi*f),E0(Epk),phi(phase),full_acceleration(true)
{}

inline void RFAcceleratingField::SetFrequency (double f)
{
    w=twoPi*f;
}

inline double RFAcceleratingField::GetFrequency () const
{
    return w/twoPi;
}

inline void RFAcceleratingField::SetPhase (double phase)
{
    phi=phase;
}

inline double RFAcceleratingField::GetPhase () const
{
    return phi;
}

inline void RFAcceleratingField::SetAmplitude (double Epk)
{
    E0=Epk;
}

inline double RFAcceleratingField::GetAmplitude () const
{
    return E0;
}

inline double RFAcceleratingField::GetWavelength () const
{
    using PhysicalConstants::SpeedOfLight;
    return twoPi*SpeedOfLight/w;
}

inline double RFAcceleratingField::GetK () const
{
    using PhysicalConstants::SpeedOfLight;
    return w/SpeedOfLight;
}

inline void RFAcceleratingField::SetFullAcceleration(bool full_acceln)
{
    full_acceleration = full_acceln;
}

inline bool RFAcceleratingField::FullAcceleration() const
{
    return full_acceleration;
}

#endif
