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

#ifndef BzField_h
#define BzField_h 1

#include "merlin_config.h"
// EMField
#include "AcceleratorModel/EMField.h"

//	Represents a constant magnetic field in along the local
//	z-axis (a solenoidal field.)

class BzField : public EMField
{
public:

    explicit BzField (double B);

    double GetStrength () const;

    //	Returns the magnetic field at the point x and time t.
    virtual Vector3D GetBFieldAt (const Point3D& x, double t = 0) const;

    //	Returns the electric field at the point x and time t
    virtual Vector3D GetEFieldAt (const Point3D& x, double t = 0) const;

    //	Sets the strength of the field in Tesla.
    void SetStrength (double B);

private:

    double Bz;
};

inline BzField::BzField (double B)
        : Bz(B)
{}

inline double BzField::GetStrength () const
{
    return Bz;
}

inline void BzField::SetStrength (double B)
{
    Bz=B;
}

#endif
