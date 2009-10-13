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

#ifndef RectMultipole_h
#define RectMultipole_h 1

#include "merlin_config.h"
// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"
// RectMultipoleField
#include "AcceleratorModel/StdComponent/RectMultipoleField.h"
// RectangularGeometry
#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"

class ComponentTracker;

//	Abstract base class for all standard rectangular
//	rectangular multipole magnets that are typically found
//	in accelerator systems.

class RectMultipole : public RectMultipoleField
{
public:

    //	Sets the primary field strength for this multipole.
    //	Units are Tesla/meter^n, where n is the primary pole
    //	type.
    void SetFieldStrength (double b);

    //	Returns the primary field strength for this multipole.
    //	Units are Tesla/meter^n, where n is the primary pole
    //	type.
    double GetFieldStrength () const;

    //	Returns the primary pole number for this multipole.
    int GetPrimaryPoleNo () const;

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 ();

    static const int ID;

protected:

    //	Constructor taking the id and the length of the magnet,
    //	and the definition of a single multipole component. b is
    //	the field (in Tesla) at the specified radius r0. if
    //	skew==true, then a skew multipole is constructed
    //	(default=false).
    RectMultipole (const string& id, double length, int npole, double b, double r0, bool skew = false);

    //	Constructor taking the id and the length of the magnet,
    //	and the definition of a single multipole component.
    //	Here, b is the n-th field derivative  (in Tesla/meter^n,
    //	where n is the pole number), evaluated at r=0.  if
    //	skew==true, then a skew multipole is constructed
    //	(default=false).
    RectMultipole (const string& id, double len, int np, double b, bool skew = false);

private:
    //	The primary pole for this multipole.
    int np;
};

inline void RectMultipole::SetFieldStrength (double b)
{
    GetField().SetFieldScale(b);
}

inline double RectMultipole::GetFieldStrength () const
{
    return GetField().GetFieldScale();
}

inline int RectMultipole::GetPrimaryPoleNo () const
{
    return np;
}

#endif
