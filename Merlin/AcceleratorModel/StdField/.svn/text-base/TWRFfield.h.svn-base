/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TWRFfield_h
#define TWRFfield_h 1

#include "merlin_config.h"
#include <cmath>
#include "AcceleratorModel/StdField/RFAcceleratingField.h"
#include "NumericalUtils/PhysicalConstants.h"

//	A RF travelling wave  field as found in typical
//	accelerating structures. The field (0,0,Ez) is defined as
//
//	        Ez(z,t)=E0 cos(k*z-w*t+phi).

class TWRFfield : public RFAcceleratingField
{
public:

    //	Constructor taking the frequency, peak electric field
    //	and the phase of the RF.
    TWRFfield (double f, double Epk, double phase = 0);

    //	Returns the magnetic field at the point x and time t.
    virtual Vector3D GetBFieldAt (const Point3D& x, double t = 0) const;

    //	Returns the electric field at the point x and time t
    virtual Vector3D GetEFieldAt (const Point3D& x, double t = 0) const;

    //	Returns the force due to this field on a particle of
    //	charge q with position x and velocity v at time t.
    virtual Vector3D GetForceAt (const Point3D& x, const Vector3D& v, double q, double t = 0) const;

    //	Calculate the Ez component.
    virtual double Ez (double z, double t) const;
};

inline TWRFfield::TWRFfield (double f, double Epk, double phase)
        : RFAcceleratingField(Epk,f,phase)
{}

inline double TWRFfield::Ez (double z, double t) const
{
    using PhysicalConstants::SpeedOfLight;
    return E0*cos(w*(z/SpeedOfLight-t)+phi);
}

#endif
