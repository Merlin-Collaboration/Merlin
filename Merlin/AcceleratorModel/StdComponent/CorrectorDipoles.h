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

#ifndef CorrectorDipoles_h
#define CorrectorDipoles_h 1

#include "merlin_config.h"
// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"

class ComponentTracker;

//	A horizontal corrector dipole.

class XCor : public RectMultipole
{
public:

    //	Constructor taking the identifier for the corrector, the
    //	length and the field (in Tesla).
    XCor (const string& id, double len, double B = 0);

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    static const int ID;
};

//	A vertical corrector dipole.

class YCor : public RectMultipole
{
public:

    //	Constructor taking the identifier for the corrector, the
    //	length and the field (in Tesla).
    YCor (const string& id, double len, double B = 0);

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    static const int ID;
};

inline XCor::XCor (const string& id, double len, double B)
        : RectMultipole(id,len,0,B)
{}

inline YCor::YCor (const string& id, double len, double B)
        : RectMultipole(id,len,0,B,true)
{}

#endif
