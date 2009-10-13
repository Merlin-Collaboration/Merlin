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
// $Revision: 1.6 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef AcceleratorComponent_h
#define AcceleratorComponent_h 1

#include "merlin_config.h"
#include <string>

// EMField
#include "AcceleratorModel/EMField.h"
// AcceleratorGeometry
#include "AcceleratorModel/AcceleratorGeometry.h"
// Aperture
#include "AcceleratorModel/Aperture.h"
// ModelElement
#include "AcceleratorModel/ModelElement.h"
// WakePotentials
#include "AcceleratorModel/WakePotentials.h"

class ComponentTracker;

using std::string;

//	An AcceleratorComponent represents any component which
//	can be placed in an accelerator lattice. Typically, a
//	lattice model is constructed as an ordered sequence of
//	AcceleratorComponents. AcceleratorComponents can be
//	associated with a geometry (an AcceleratorGeometry
//	object) and a em field (an EMFieldRegion object), which
//	when taken together uniquely define the field properties
//	for the component. Component "tracking" is supported via
//	the funtion PrepareTracker(Tracker&), which sets up a
//	Tracker object to track the component.

class AcceleratorComponent : public ModelElement
{
public:
    virtual ~AcceleratorComponent ();

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the geometry length of the component.
    double GetLength () const;

    //	Returns a pointer to the this components geometry.
    //	Returns NULL if no geometry is associated with this
    //	component.
    const AcceleratorGeometry* GetGeometry () const;

    //	Returns a pointer to this components field. A NULL
    //	pointer is returned if the component has no field.
    const EMField* GetEMField () const;

    const Aperture* GetAperture () const;

    void SetAperture (Aperture* ap);

    //	Returns the wake potentials associated with this cavity.
    WakePotentials* GetWakePotentials () const;

    //	Sets the wake potentials associated with this cavity.
    void SetWakePotentials (WakePotentials* wp);

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Rotates the component 180 degrees about its local Y axis.
    virtual void RotateY180 () = 0;

	//  Set/Get the uniques beamline index for this frame
	void SetBeamlineIndex(size_t n);
	size_t GetBeamlineIndex() const;
	void AppendBeamlineIndecies(std::vector<size_t>&) const;

    //	Returns the total number of distinct component types in
    //	the system.
    static int TotalComponentNumber ();

    // Data Members for Class Attributes

    //	Unique index for an Accelerator component.
    static const int ID;

protected:

    //	Protected constructors used by derived classes.
    explicit AcceleratorComponent (const string& aName = string());

    AcceleratorComponent (const string& aName, AcceleratorGeometry* aGeom, EMField* aField);
	

    //	Used by derived classes to generate a unique index. All
    //	derived classes should have a static member ID of type
    //	IndexType which should be initialised as follows
    //
    //	IndexType component::ID = UniqueIndex();
    static int UniqueIndex ();

    EMField* itsField;
    AcceleratorGeometry* itsGeometry;
    Aperture* itsAperture;
    WakePotentials* itsWakes;

	// beamline index associated with this component
	size_t blI; 
};


inline AcceleratorComponent::AcceleratorComponent (const string& aName)
        : ModelElement(aName),itsField(0),itsGeometry(0),itsAperture(0),itsWakes(0)
{}

inline AcceleratorComponent::AcceleratorComponent (const string& aName, AcceleratorGeometry* aGeom, EMField* aField)
        : ModelElement(aName),itsField(aField),itsGeometry(aGeom),itsAperture(0),itsWakes(0)
{}

inline const AcceleratorGeometry* AcceleratorComponent::GetGeometry () const
{
    return itsGeometry;
}

inline const EMField* AcceleratorComponent::GetEMField () const
{
    return itsField;
}

inline const Aperture* AcceleratorComponent::GetAperture () const
{
    return itsAperture;
}

inline void AcceleratorComponent::SetAperture (Aperture* ap)
{
    itsAperture = ap;
}

inline WakePotentials* AcceleratorComponent::GetWakePotentials () const
{
    return itsWakes;
}

inline void AcceleratorComponent::SetWakePotentials (WakePotentials* wp)
{
    itsWakes=wp;
}

inline int AcceleratorComponent::TotalComponentNumber ()
{
    return 0;
}

inline void AcceleratorComponent::SetBeamlineIndex(size_t n)
{
	blI = n;
}

inline size_t AcceleratorComponent::GetBeamlineIndex() const
{
	return blI;
}

inline void AcceleratorComponent::AppendBeamlineIndecies(std::vector<size_t>& ivec) const
{
	ivec.push_back(blI);
}

#endif
