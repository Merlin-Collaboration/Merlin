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

#ifndef Drift_h
#define Drift_h 1

#include "merlin_config.h"
// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// SimpleDrift
#include "AcceleratorModel/StdComponent/SimpleDrift.h"
// RectangularGeometry
#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"

class ComponentTracker;

//	A simple drift section.

class Drift : public SimpleDrift
{
public:

    explicit Drift (const string& id, double len = 0);

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

    //	Unique index for an Accelerator component.
    static const int ID;
};


#endif
