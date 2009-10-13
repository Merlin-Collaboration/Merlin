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

#ifndef SWRFStructure_h
#define SWRFStructure_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdComponent/RFStructure.h"

class ComponentTracker;

//	A standing wave accelerating structure.

class SWRFStructure : public RFStructure
{
public:
    //	Constructor taking the name of the cavity, the number of
    //	cells, the frequency (f)  in MHz, and peak electric
    //	field (E0) in MV/m and the phase.
    SWRFStructure (const string& id, int ncells, double f, double E0, double phi = 0);

    //	Returns the type string for this component.
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

    // Data Members for Class Attributes

    static const int ID;

protected:
private:
	// Private copy constuctor
	SWRFStructure (const SWRFStructure&);
};

// Class SWRFStructure



#endif
