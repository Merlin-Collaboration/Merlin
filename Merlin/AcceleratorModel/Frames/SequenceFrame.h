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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef SequenceFrame_h
#define SequenceFrame_h 1

#include "merlin_config.h"
#include <list>
#include <vector>

// LatticeFrame
#include "AcceleratorModel/Frames/LatticeFrame.h"

class SequenceGeometry;

#define DEF_SEQUENCE_LABEL "<UNNAMED>"

//	A  LatticeFrame object which is composed of a contiguous
//	sequence of other (sub-)LatticeFrame objects.

class SequenceFrame : public LatticeFrame
{
public:

    //	Enumeration constants which define location of origin
    //	for this frame.
    typedef enum {originAtEntrance,originAtCenter,originAtExit} Origin;
    typedef std::list<LatticeFrame*> FrameList;

    explicit SequenceFrame (const string& id = DEF_SEQUENCE_LABEL, Origin originLoc = originAtCenter);
    //	Copy constructor
    SequenceFrame (const SequenceFrame& rhs);

    virtual ~SequenceFrame ();

    // Propagate FrameTraverser object (iteration)
    virtual void Traverse(FrameTraverser& ft);

    //	Sequence construction: append aFrame to the sequence.
    void AppendFrame (LatticeFrame& aFrame);

    //	Function called after construction of the Accelerator
    //	Model is complete. Allows the nested frame hierachy to
    //	perform certain state checks and updates, which are only
    //	possible once the entire model is complete.
    virtual void ConsolidateConstruction ();

    //	Causes any cached state to be invalidated. The cached
    //	state should be re-calculated if and when required.
    virtual void Invalidate () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Replace subFrame with newSubFrame. Returns true if
    //	successful (i.e. subFrame is a sub-frame of this Lattice
    //	Frame).
    virtual bool ReplaceSubFrame (LatticeFrame* subFrame, LatticeFrame* newSubFrame);

    //	Return the type string for the element.
    virtual const string& GetType () const;

	// constucting bealine index list
	void AppendBeamlineIndecies(std::vector<size_t>&) const;

protected:


private:

    FrameList subFrames;
    SequenceGeometry* itsSeqGeom;

    //	Copies the subframes from frames.
    void CopySubFrames (const list<LatticeFrame*>& frames);

    //	Returns true if aSubFrame is the first (entrance)
    //	sub-frame of this frame.
    virtual bool IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const;

    // don't allow copy assignment
    const SequenceFrame & operator=(const SequenceFrame &right);
};

#endif
