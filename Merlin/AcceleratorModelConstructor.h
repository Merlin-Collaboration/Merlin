/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorModelConstructor_h
#define AcceleratorModelConstructor_h 1

#include "merlin_config.h"
#include <stack>
#include "TComponentFrame.h"
#include "AcceleratorModel.h"
#include "utils.h"

class SequenceFrame;

/**
 *	Responsible for constructing an AcceleratorModel object.
 *	The nested frame hierarchy is constructed internally
 *	using a frame stack: each call to NewFrame(aFrame)
 *	pushes the current frame onto the stack, and makes the
 *	new frame the current frame. A subsequent call to End
 *	Frame(aFrame) causes the current frame (aFrame) to be
 *	popped from the stack. At any time between calls to New
 *	Frame and EndFrame can AppendComponent be called.
 */
class AcceleratorModelConstructor
{
public:

	typedef std::stack<SequenceFrame*> FrameStack;

	AcceleratorModelConstructor();
	~AcceleratorModelConstructor();

	/**
	 *	Initialises a new AcceleratorModel. Must be called
	 *	before any subsequent constructor calls.
	 */
	void NewModel();

	/**
	 *	Ends the current model construction and returns the
	 *	model. The model must be in a valid (complete) state.
	 */
	AcceleratorModel* GetModel();

	/**
	 *	Begin construction of a new LatticeFrame.
	 */
	void NewFrame(SequenceFrame* aFrame);

	/**
	 *	Finish construction of the specified frame.
	 */
	void EndFrame();

	/*
	 *	Append the specified component a distance d meters
	 *	downstream of the last component.
	 */
	template<class T> T* AppendComponent(T& acc, double d = 0)
	{
		if(!fequal(d, 0.0))
		{
			AppendDrift(d);
		}
		AppendComponentFrame(new TComponentFrame<T>(acc));
		return &acc;
	}

	template<class T> T* AppendComponent(T* acc, double d = 0)
	{
		return AppendComponent(*acc, d);
	}

	/**
	 * Append an arbitrary ComponentFrame
	 */
	void AppendComponentFrame(ComponentFrame* cf);

	/*
	 * Append an entire SequenceFrame (or derivative) to the current model.
	 * This function allows complex structures (e.g. girders) which have
	 * been externally constructed to placed in the current model
	 */
	void AppendFrame(SequenceFrame* aFrame);

	/**
	 *	Returns a const reference to the current frame.
	 *	@return The current frame reference
	 */
	SequenceFrame& GetCurrentFrame() const
	{
		return *(frameStack.top());
	}

	/**
	 *	Returns the depth of the current frame. 0 refers to the
	 *	global frame (top level).
	 *	@return The current frame depth
	 */
	int GetCurrentFrameDepth() const
	{
		return frameStack.size() - 1;
	}

	/**
	 *	Appends a simple drift to the current model.
	 */
	void AppendDrift(double d);

	/**
	 *	Adds a ModelElement to the Accelerator Model.
	 */
	void AddModelElement(ModelElement* element);

	/**
	 *	Prints a table to os containing statistics on the type
	 *	and number of ModelElement current contained in the
	 *	model.
	 */
	void ReportStatistics(std::ostream& os) const;

private:

	AcceleratorModel* currentModel;
	FrameStack frameStack;
};

#endif
