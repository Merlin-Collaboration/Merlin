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

#ifndef Marker_h
#define Marker_h 1

#include "merlin_config.h"
// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"

class ComponentTracker;

//	Special Marker component. This does not represent a
//	physical component, and is provided for uses wishing to
//	"mark" specific locations in the lattice. It has no
//	field or geometry associated with it.

class Marker : public AcceleratorComponent
{
public:

    //	Constructor
    explicit Marker (const std::string& id);

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

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();

    //	Unique index for an Accelerator component.
    static const int ID;
};

inline Marker::Marker (const std::string& id)
        : AcceleratorComponent(id,0,0)
{}

#endif
