/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef Measurement_h
#define Measurement_h 1

//	POD representing a physically measured quantity, which
//	has a value and an associated error.

struct Measurement
{
    Measurement (double v, double err);
    Measurement ();

    double value;
    double error;
};

// Class Measurement

inline Measurement::Measurement (double v, double err)
        : value(v),error(err)
{}

inline Measurement::Measurement ()
{}

#endif
