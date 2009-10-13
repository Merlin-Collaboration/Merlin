//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\AcceleratorModel\StdComponent\Solenoid.h
 * last modified 10/12/01 16:41:41
 */

/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 2001 by The Merlin Collaboration.
 * - ALL RIGHTS RESERVED - 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */


#ifndef Solenoid_h
#define Solenoid_h 1

#include "merlin_config.h"


// SimpleSolenoid
#include "AcceleratorModel/StdComponent/SimpleSolenoid.h"


//	A simple solenoid with field Bz.



class Solenoid : public SimpleSolenoid
{
public:
    Solenoid (const std::string& id, double len, double Bz);


    //	Returns the value of the field in Tesla.
    double GetBz () const;

    //	Sets the value of the field in Tesla.
    void SetBz (double B);

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();

    //	Return the type string for the element.
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

	// The followind field access function added for
	// compatability with other magnets
	void SetFieldStrength(double b) { SetBz(b); }
	double GetFieldStrength() const { return GetBz(); }

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};

// Class Solenoid


inline double Solenoid::GetBz () const
{
    return GetField().GetStrength();
}

inline void Solenoid::SetBz (double B)
{
    GetField().SetStrength(B);
}



#endif
