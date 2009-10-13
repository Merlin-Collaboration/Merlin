/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:13:00 $
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ComponentFrame_h
#define ComponentFrame_h 1

#include "merlin_config.h"
// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"
// LatticeFrame
#include "AcceleratorModel/Frames/LatticeFrame.h"

//	An atomic LatticeFrame object which represents a single
//	instance of an AcceleratorComponent in a beamline.

class ComponentFrame : public LatticeFrame
{
public:

    //	Constructor taking the associated AcceleratorComponent
    //	object.
    explicit ComponentFrame (AcceleratorComponent& ac, const string& id = "");

    //	Copy constructor.
    ComponentFrame (const ComponentFrame& rhs);

    virtual ~ComponentFrame ();

    AcceleratorComponent& GetComponent ();
    const AcceleratorComponent& GetComponent () const;

    // Returns true if this is an empty frame
    bool IsComponent() const {
        return theComponent!=0;
    }

    // Returns the (design) geometry patches for this component
    // frame.
    virtual const Transform3D* GetEntranceGeometryPatch() const {
        return 0;
    }
    virtual const Transform3D* GetExitGeometryPatch() const {
        return 0;
    }

    //	Causes any cached state to be invalidated. The cached
    //	state should be re-calculated if and when required.
    virtual void Invalidate () const;

    //	Return the name of the element. Returns the name of the
    //	AcceleratorComponent if the label for this frame has not
    //	been explicitely set.
    virtual const string& GetName () const;

    //	Returns the type of the referenced AcceleratorComponent.
    virtual const string& GetType () const;

    //	Returns a copy of this ComponentFrame. Note that only
    //	the reference to the AcceleratorComponent is copied, not
    //	the AcceleratorComponent itself.
    virtual ModelElement* Copy () const;

	//  Set/Get the uniques beamline index for this frame
	void SetBeamlineIndex(size_t n);
	size_t GetBeamlineIndex() const;
	void AppendBeamlineIndecies(std::vector<size_t>&) const;

protected:

    ComponentFrame(AcceleratorComponent* ac, const string& id = "");

    AcceleratorComponent* theComponent;

    //	Should  never be called.
    virtual bool IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const;


private:
	size_t blI; // beamline index
};

inline ComponentFrame::ComponentFrame (AcceleratorComponent& ac, const string& id)
        : LatticeFrame(id.empty()?ac.GetQualifiedName():id),theComponent(&ac)
{
    SetGeometry(theComponent->GetGeometry());
}

inline ComponentFrame::ComponentFrame (const ComponentFrame& rhs)
: LatticeFrame(rhs),theComponent(rhs.theComponent)
{}

inline ComponentFrame::ComponentFrame(AcceleratorComponent* ac, const string& id)
: LatticeFrame(id),theComponent(ac)
{}

inline AcceleratorComponent& ComponentFrame::GetComponent ()
{
    if(theComponent==0) throw "bad_component";
    return *theComponent;
}

inline const AcceleratorComponent& ComponentFrame::GetComponent () const
{
    if(theComponent==0) throw "bad_component";
    return *theComponent;
}

inline const string& ComponentFrame::GetName () const
{
    const string& id = LatticeFrame::GetName();
    return (id.length()!=0)? id : (theComponent->GetName());
}

inline void ComponentFrame::SetBeamlineIndex(size_t n)
{
	blI = n;
	if(theComponent)
		theComponent->SetBeamlineIndex(n);
}

inline size_t ComponentFrame::GetBeamlineIndex() const
{
	return blI;
}

inline void ComponentFrame::AppendBeamlineIndecies(std::vector<size_t>& ivec) const
{
	ivec.push_back(blI);
}

#endif
