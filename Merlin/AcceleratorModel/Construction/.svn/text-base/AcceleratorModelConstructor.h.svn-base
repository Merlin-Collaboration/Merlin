/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:33:05 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef AcceleratorModelConstructor_h
#define AcceleratorModelConstructor_h 1

#include "merlin_config.h"
#include <stack>
// TComponentFrame
#include "AcceleratorModel/Frames/TComponentFrame.h"
// AcceleratorModel
#include "AcceleratorModel/AcceleratorModel.h"

class SequenceFrame;

//	Responsible for constructing an AcceleratorModel object.
//	The nested frame hierachy is constructed internally
//	using a frame stack: each call to NewFrame(aFrame)
//	pushes the current frame onto the stack, and makes the
//	new frame the current frame. A subsequent call to End
//	Frame(aFrame) causes the current frame (aFrame) to be
//	popped from the stack. At any time between calls to New
//	Frame and EndFrame can AppendComponent be called.

class AcceleratorModelConstructor
{
public:

    typedef std::stack< SequenceFrame* > FrameStack;

    AcceleratorModelConstructor ();
    ~AcceleratorModelConstructor ();

    //	Initialises a new AcceleratorModel. Must be called
    //	before any subsequent constructor calls.
    void NewModel ();

    //	Ends the current model construction and returns the
    //	model. The model must be in a valid (complete) state.
    AcceleratorModel* GetModel ();

    //	Begin construction of a new LatticeFrame.
    void NewFrame (SequenceFrame* aFrame);

    //	Finish construction of the specified frame.
    void EndFrame ();

    //	Append the specified component a distance d meters
    //	downstream of the last component.
    template<class T> T* AppendComponent (T& acc, double d = 0)
    {
        if(d!=0)
            AppendDrift(d);
        AppendComponentFrame(new TComponentFrame<T>(acc));
        return &acc;
    }

    template<class T> T* AppendComponent(T* acc, double d=0)
    {
        return AppendComponent(*acc,d);
    }

    // Append an arbitrary ComponentFrame
    void AppendComponentFrame (ComponentFrame* cf);

    // Append an entire SequenceFrame (or derivative) to the current model.
    // This function allows complex structures (eg. girders) which have
    // been externally constructed to placed in the current model
    void AppendFrame(SequenceFrame* aFrame);

    //	Returns a const reference to the current frame.
    SequenceFrame& GetCurrentFrame () const
    {
        return *(frameStack.top());
    }

    //	Returns the depth of the current frame. 0 refers to the
    //	global frame (top level).
    int GetCurrentFrameDepth () const
    {
        return frameStack.size()-1;
    }

    //	Appends a simple drift to the current model.
    void AppendDrift (double d);

    //	Adds a ModelElement to the Accelerator Model.
    void AddModelElement (ModelElement* element);

    //	Prints a table to os  containing statistics on the type
    //	and number of ModelElement current contained in the
    //	model.
    void ReportStatistics (std::ostream& os) const;

private:

    AcceleratorModel* currentModel;
    FrameStack frameStack;
};


#endif
