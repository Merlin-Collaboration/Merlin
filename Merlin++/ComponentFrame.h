/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ComponentFrame_h
#define ComponentFrame_h 1

#include "merlin_config.h"
#include "AcceleratorComponent.h"
#include "LatticeFrame.h"

/**
 *	An atomic LatticeFrame object which represents a single
 *	instance of an AcceleratorComponent in a beamline.
 */

class ComponentFrame: public LatticeFrame
{
public:

	/**
	 *	Constructor taking the associated AcceleratorComponent
	 *	object.
	 */
	explicit ComponentFrame(AcceleratorComponent& ac, const string& id = "");

	/**
	 *	Copy constructor.
	 */
	ComponentFrame(const ComponentFrame& rhs);

	virtual ~ComponentFrame();

	AcceleratorComponent& GetComponent();
	const AcceleratorComponent& GetComponent() const;

	/**
	 * Returns true if this is an empty frame
	 */
	bool IsComponent() const
	{
		return theComponent != nullptr;
	}

	/**
	 * Returns the (design) geometry patches for this component
	 * frame.
	 */
	virtual const Transform3D* GetEntranceGeometryPatch() const
	{
		return nullptr;
	}
	virtual const Transform3D* GetExitGeometryPatch() const
	{
		return nullptr;
	}

	/**
	 *	Causes any cached state to be invalidated. The cached
	 *	state should be re-calculated if and when required.
	 */
	virtual void Invalidate() const;

	/**
	 *	Return the name of the element. Returns the name of the
	 *	AcceleratorComponent if the label for this frame has not
	 *	been explicitly set.
	 */
	virtual const string& GetName() const;

	/**
	 *	Returns the type of the referenced AcceleratorComponent.
	 *	@return Referenced AcceleratorComponent type
	 */
	virtual const string& GetType() const;

	/**
	 *	Returns a copy of this ComponentFrame. Note that only
	 *	the reference to the AcceleratorComponent is copied, not
	 *	the AcceleratorComponent itself.
	 *	@return Copy of the ComponentFrame
	 */
	virtual ModelElement* Copy() const;

	/**
	 * Set the unique beamline index for this frame
	 */
	void SetBeamlineIndex(size_t n);

	/**
	 * Get the unique beamline index for this frame
	 */
	size_t GetBeamlineIndex() const;
	void AppendBeamlineIndexes(std::vector<size_t>&) const;

protected:

	ComponentFrame(AcceleratorComponent* ac, const string& id = "");

	AcceleratorComponent* theComponent;

	//	Should  never be called.
	virtual bool IsBoundaryPlane(BoundaryPlane p, const LatticeFrame* aSubFrame) const;

private:
	size_t blI; /// beamline index

	//Just need to disable the assignment operator - copy constructor is defined
	ComponentFrame& operator=(const ComponentFrame& frame);

};

inline ComponentFrame::ComponentFrame(AcceleratorComponent& ac, const string& id) :
	LatticeFrame(id.empty() ? ac.GetQualifiedName() : id), theComponent(&ac), blI(0)
{
	SetGeometry(theComponent->GetGeometry());
}

inline ComponentFrame::ComponentFrame(const ComponentFrame& rhs) :
	LatticeFrame(rhs), theComponent(rhs.theComponent), blI(0)
{
}

inline ComponentFrame::ComponentFrame(AcceleratorComponent* ac, const string& id) :
	LatticeFrame(id), theComponent(ac), blI(0)
{
}

inline AcceleratorComponent& ComponentFrame::GetComponent()
{
	if(theComponent == nullptr)
	{
		throw "bad_component";
	}
	return *theComponent;
}

inline const AcceleratorComponent& ComponentFrame::GetComponent() const
{
	if(theComponent == nullptr)
	{
		throw "bad_component";
	}
	return *theComponent;
}

inline const string& ComponentFrame::GetName() const
{
	const string& id = LatticeFrame::GetName();
	return (id.length() != 0) ? id : (theComponent->GetName());
}

inline void ComponentFrame::SetBeamlineIndex(size_t n)
{
	blI = n;
	if(theComponent)
	{
		theComponent->SetBeamlineIndex(n);
	}
}

inline size_t ComponentFrame::GetBeamlineIndex() const
{
	return blI;
}

inline void ComponentFrame::AppendBeamlineIndexes(std::vector<size_t>& ivec) const
{
	ivec.push_back(blI);
}

#endif
