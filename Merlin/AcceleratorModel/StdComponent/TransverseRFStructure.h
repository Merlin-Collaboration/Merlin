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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TransverseRFStructure_h
#define TransverseRFStructure_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdComponent/RFStructure.h"
#include "AcceleratorModel/StdField/TransRFfield.h"

class ComponentTracker;

//	A travelling wave transverse deflecting structure.

class TransverseRFStructure : public RFStructure
{
public:

    //	Constructor taking the label for this cavity (id), the
    //	cavity length (len) in meters, the frequency (f) in MHz,
    //	and gradient (Epk) in MV/m and the phase (phi) in
    //	radians.
    TransverseRFStructure (const string& id, double len, double f, double Epk, double phi=0, double theta=0);

    //  Orientation
    double GetFieldOrientation() const;
    void SetFieldOrientation(double);

    //	Returns the type string "TransverseRFStructure".
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

protected:
private:
	// Private copy constructor
	TransverseRFStructure(const TransverseRFStructure&);
};

inline double TransverseRFStructure::GetFieldOrientation() const
{
    return static_cast<TransverseRFfield*>(itsField)->GetFieldOrientation();
}

inline void TransverseRFStructure::SetFieldOrientation(double t)
{
    static_cast<TransverseRFfield*>(itsField)->SetFieldOrientation(t);
}

#endif
