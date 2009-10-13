//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\AcceleratorModel\StdComponent\StandardMultipoles.h
 * last modified 16/05/02 11:10:34
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


#ifndef StandardMultipoles_h
#define StandardMultipoles_h 1

#include "merlin_config.h"


// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"

class ComponentTracker;

// type strings for template parameter instantiation


//	A standard quadrupole magnet.




class Quadrupole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the quadrupole gradient. Units are Tesla/meter.
    Quadrupole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    Quadrupole (const string& id, double len, double B, double r0);


    //	Initialise aTracker with *this.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Returns the unique index for a Quadrupole
    virtual int GetIndex () const;

    //	Returns the type string "Quadrupole"
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};

//	A standard sextupole magnet.




class Sextupole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the sextupole gradient. Units are Tesla/meter^2.
    Sextupole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    Sextupole (const string& id, double len, double B, double r0);


    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};

//	 A standard skew-quadrupole magnet.




class SkewQuadrupole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the quadrupole gradient. Units are Tesla/meter.
    SkewQuadrupole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    SkewQuadrupole (const string& id, double len, double B, double r0);


    //	Initialise aTracker to track *this.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Return the index for a SkewQuadrupole.
    virtual int GetIndex () const;

    //	Returns the type string "SkewQuadrupole".
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};

//	A standard Octupole magnet.




class Octupole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the octupole gradient. Units are Tesla/meter^3.
    Octupole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    Octupole (const string& id, double len, double B, double r0);


    //	Initialise aTracker to track *this.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Returns the index for an Octupole.
    virtual int GetIndex () const;

    //	Returns the type string "Octupole".
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};




class SkewSextupole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the sextupole gradient. Units are Tesla/meter^2.
    SkewSextupole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    SkewSextupole (const string& id, double len, double B, double r0);


    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
private:
};

//	A standard (normal) decapole magnet.



class Decapole : public RectMultipole
{
public:
    //	Constructor. The primary field component is specified as
    //	a the octupole gradient. Units are Tesla/meter^4.
    Decapole (const string& id, double len, double dnB);

    //	Constructor. The primary field component is specified as
    //	a field B (in Tesla) at a speciific radius r0 (in meter).
    Decapole (const string& id, double len, double B, double r0);


    //	Initialise aTracker to track *this.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Returns the index for an Octupole.
    virtual int GetIndex () const;

    //	Returns the type string "Octupole".
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    static const int ID;

protected:
private:
private:
};

// Class Quadrupole

// Class Sextupole

// Class SkewQuadrupole

// Class Octupole

// Class SkewSextupole

// Class Decapole



#endif
