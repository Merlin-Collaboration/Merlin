/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorComponent_h
#define AcceleratorComponent_h 1

#include "merlin_config.h"
#include <string>

#include "EMField.h"
#include "AcceleratorGeometry.h"
#include "Aperture.h"
#include "ModelElement.h"
#include "WakePotentials.h"

class ComponentTracker;

using std::string;

/**
 * An AcceleratorComponent represents any component which
 * can be placed in an accelerator lattice. Typically, a
 * lattice model is constructed as an ordered sequence of
 * AcceleratorComponents. AcceleratorComponents can be
 * associated with a geometry (an AcceleratorGeometry
 * object) and a em field (an EMField Region object), which
 * when taken together uniquely define the field properties
 * for the component. Component "tracking" is supported via
 * the function PrepareTracker(Tracker&), which sets up a
 * Tracker object to track the component.
 */
class AcceleratorComponent: public ModelElement
{
public:

	virtual ~AcceleratorComponent();

	/**
	 * Returns the unique index for this class of accelerator
	 * components.
	 * @return An integer containing the unique index for this AcceleratorComponent type.
	 */
	virtual int GetIndex() const;

	/**
	 * Returns the geometry length of the component.
	 * @return A double containing the length of this component.
	 */
	double GetLength() const;

	/**
	 * Returns a pointer to the this components geometry.
	 * Returns a nullptr if no geometry is associated with this component.
	 * @return The AcceleratorGeometry associated with this element.
	 */
	const AcceleratorGeometry* GetGeometry() const;

	/**
	 * Returns a pointer to this components field. A nullptr is returned if the component has no field.
	 * @return The EMField associated with this element.
	 */
	const EMField* GetEMField() const;

	/**
	 * Returns a pointer to the aperture of this element.
	 * @return The Aperture associated with this element.
	 */
	Aperture* GetAperture() const;

	/**
	 * Sets the aperture of this element
	 * @param[in] ap A pointer to an Aperture class to associate with this component
	 */
	void SetAperture(Aperture* ap);

	/**
	 * Returns the wake potentials associated with this element.
	 * @return The WakePotentials associated with this element.
	 */
	WakePotentials* GetWakePotentials() const;

	/**
	 * Sets the wake potentials associated with this cavity.
	 * @param[in] wp A pointer to a WakePotentials class to associate with this component
	 */
	void SetWakePotentials(WakePotentials* wp);

	/**
	 * Primary tracking interface. Prepares the specified
	 * Tracker object for tracking this component.
	 * @param[in,out] aTracker The tracker to prepare.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	/**
	 * Rotates the component 180 degrees about its local Y axis.
	 */
	virtual void RotateY180() = 0;

	/**
	 * Set the unique beamline index for this frame.
	 * @param[in] n The beamline index.
	 */
	void SetBeamlineIndex(size_t n);

	/**
	 * Get the unique beamline index for this frame.
	 * @return A size_t containing the unique beamline index.
	 */
	size_t GetBeamlineIndex() const;

	void AppendBeamlineIndexes(std::vector<size_t>&) const;

	/**
	 * Returns the total number of distinct component types in
	 * the system.
	 * @return An integer containing the total number of types of element types that exist.
	 */
	static int TotalComponentNumber();

	// Data Members for Class Attributes
	/**
	 * Unique index for an Accelerator component.
	 */
	static const int ID;

	/**
	 * Set the distance from the start of the lattice to the START of the element.
	 * @param[in] position The position of the start of the lattice element.
	 */
	void SetComponentLatticePosition(double position);

	/**
	 * Get the distance from the start of the lattice to the START of the element.
	 * @return A double containing the position of the start of the lattice element.
	 */
	double GetComponentLatticePosition() const;

	/**
	 * Collimator ID for FLUKA output AV+HR 09.11.15
	 * N.B. These values can only be set or got in Collimator.
	 * @param[in] n The value of the Collimator ID to set.
	 */
	void SetCollID(int n)
	{
		Coll_ID = n;
	}

	/**
	 * Get the Collimator ID.
	 * @return An integer containing the Collimator ID.
	 */
	int GetCollID() const
	{
		return Coll_ID;
	}

protected:

	/**
	 * Protected constructors used by derived classes.
	 * @param[in] aName The name of the AcceleratorComponent
	 */
	explicit AcceleratorComponent(const string& aName = string());

	/**
	 * Protected constructors used by derived classes.
	 * @param[in] aName The name of the AcceleratorComponent
	 * @param[in] aGeom The AcceleratorGeometry of the AcceleratorComponent
	 * @param[in] aField The EMField of the AcceleratorComponent
	 */
	AcceleratorComponent(const string& aName, AcceleratorGeometry* aGeom, EMField* aField);

	/**
	 * Used by derived classes to generate a unique index. All
	 * derived classes should have a static member ID of type
	 * IndexType which should be initialised as follows:
	 *
	 * IndexType component::ID = UniqueIndex();
	 */
	static int UniqueIndex();

	/**
	 * A pointer to the Electromagnetic field for this component
	 */
	EMField* itsField;

	/**
	 * A pointer to the geometry for this component
	 */
	AcceleratorGeometry* itsGeometry;

	/**
	 * A pointer to the Aperture for this component
	 */
	Aperture* itsAperture;

	/**
	 * A pointer to the Wake potentials for this component
	 */
	WakePotentials* itsWakes;

	/**
	 * The position of this element relative to a user defined location
	 * (the start of the user generated lattice).
	 */
	double position;

	/**
	 * beamline index associated with this component
	 */
	size_t blI;

	/**
	 * Collimator ID for FLUKA output AV+HR 09.11.15
	 * N.B. These values can only be set or got in Collimator
	 */
	int Coll_ID;

};

inline AcceleratorComponent::AcceleratorComponent(const string& aName) :
	ModelElement(aName), itsField(nullptr), itsGeometry(nullptr), itsAperture(nullptr), itsWakes(nullptr), position(0),
	blI(0)
{
}

inline AcceleratorComponent::AcceleratorComponent(const string& aName, AcceleratorGeometry* aGeom, EMField* aField) :
	ModelElement(aName), itsField(aField), itsGeometry(aGeom), itsAperture(nullptr), itsWakes(nullptr), position(0),
	blI(0)
{
}

inline const AcceleratorGeometry* AcceleratorComponent::GetGeometry() const
{
	return itsGeometry;
}

inline const EMField* AcceleratorComponent::GetEMField() const
{
	return itsField;
}

inline Aperture* AcceleratorComponent::GetAperture() const
{
	return itsAperture;
}

inline void AcceleratorComponent::SetAperture(Aperture* ap)
{
	itsAperture = ap;
}

inline WakePotentials* AcceleratorComponent::GetWakePotentials() const
{
	return itsWakes;
}

inline void AcceleratorComponent::SetWakePotentials(WakePotentials* wp)
{
	itsWakes = wp;

}

inline int AcceleratorComponent::TotalComponentNumber()
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

inline void AcceleratorComponent::AppendBeamlineIndexes(std::vector<size_t>& ivec) const
{
	ivec.push_back(blI);
}

inline void AcceleratorComponent::SetComponentLatticePosition(double distance)
{
	position = distance;
}

inline double AcceleratorComponent::GetComponentLatticePosition() const
{
	return position;
}

#endif
