/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		06.10.14 Haroon Rafique
// Modified:		
// Last Edited: 05.01.15 HR
// 
/////////////////////////////////////////////////////////////////////////
#ifndef HollowElectronLens_h
#define HollowElectronLens_h 1

#include "merlin_config.h"

#include "AcceleratorModel/StdComponent/Drift.h"

class ComponentTracker;

// A HollowElectronLens provides a method of active halo control.
// It can be used as a soft scraper. Optically the element is a drift.
// See HollowELensProcess.
 
class HollowElectronLens : public Drift 
{
public:

	HollowElectronLens (const string& id, double len);

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
    
    virtual double GetRmax () const {return Rmax;}
    
    virtual double GetRmin () const {return Rmin;}
    
    virtual void SetRmax (double rmax);
    
    virtual void SetRmin (double rmin);

    //	Unique index for an Accelerator component.
    static const int ID;

private:

	double Rmin;
	double Rmax;

};


#endif
