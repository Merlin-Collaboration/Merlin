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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TWRFStructure_h
#define TWRFStructure_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdComponent/RFStructure.h"

class ComponentTracker;

//	A travelling wave accelerating structure.

class TWRFStructure : public RFStructure
{
public:

    //	Constructor taking the label for this cavity (id), the
    //	cavity length (len) in meters, the frequency (f) in MHz,
    //	and gradient (Epk) in MV/m and the phase (phi) in
    //	radians.
    TWRFStructure (const string& id, double len, double f, double Epk, double phi = 0);
    TWRFStructure (const TWRFStructure&);

    //	Returns the type string "TWRFStructure".
    virtual const string& GetType () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    static const int ID;
};


#endif
