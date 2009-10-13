/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef FrameModifier_h
#define FrameModifier_h 1

#include "merlin_config.h"
// LatticeFrame
#include "AcceleratorModel/Frames/LatticeFrame.h"

//	A utility class which adds an additional layer to the
//	frame hierachy, and thus another level of coordinate
//	transformations. FrameModifier effectively "wraps" a
//	single LatticeFrame object. It can insert and remove
//	itself in the lattice frame hierachy.

class FrameModifier : public LatticeFrame
{
public:

    //	Constructor taking the LatticeFrame object above which
    //	this is to be inserted.
    explicit FrameModifier (LatticeFrame* frame, const std::string& label ="");
    ~FrameModifier();

    //	Returns a copy of the sub-frame.
    ModelElement* Copy () const;

    //	Returns the type string of the sub-frame.
    const string& GetType () const;

    //	Causes any cached state to be invalidated. The cached
    //	state should be re-calculated if and when required.
    virtual void Invalidate () const;

    //	Function called after construction of the Accelerator
    //	Model is complete. Allows the nested frame hierachy to
    //	perform certain state checks and updates, which are only
    //	possible once the entire model is complete.
    virtual void ConsolidateConstruction ();

    // Remove the modifier frame from the hierachy
    void Remove();

private:

    LatticeFrame* subFrame;

    virtual bool IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const;
};

#endif
