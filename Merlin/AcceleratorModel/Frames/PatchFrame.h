/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/05/10 14:00:07 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef PatchFrame_h
#define PatchFrame_h 1

#include "merlin_config.h"
// ComponentFrame
#include "AcceleratorModel/Frames/ComponentFrame.h"
// PatchGeometry
#include "AcceleratorModel/StdGeometry/GeometryPatch.h"

// A special ComponentFrame which represents a pure geometry patch.
// There is no AcceleratorComponent associated with this frame.
//
// NOTE*** This his is a bit of a kludge, and there is probably a better
//         way to do this!

class PatchFrame : public ComponentFrame
{
public:

    //	Constructor  
    explicit PatchFrame (GeometryPatch* pg, const string& id = "");

    //	Copy constructor.
    PatchFrame (const PatchFrame& rhs);

    // Destructor
    virtual ~PatchFrame ();
    virtual const string& GetType () const {
        _TYPESTR(PatchFrame)
    }

    virtual const string& GetName () const {
        return LatticeFrame::GetName();
    }

    virtual ModelElement* Copy () const {
        return new PatchFrame(*this);
    }

    virtual const Transform3D* GetEntranceGeometryPatch() const {
        return itsPatch ? itsPatch->GetTransformation() : 0;
    }

    void SetGeometryPatch(GeometryPatch* gp) {
        if(itsPatch)
            delete itsPatch;
        itsPatch = gp;
        SetGeometry(itsPatch);
    }

private:

    GeometryPatch* itsPatch;
};

inline PatchFrame::PatchFrame (GeometryPatch* pg, const string& id)
        : ComponentFrame(0,id), itsPatch(pg)
{
    SetGeometry(itsPatch);
}

inline PatchFrame::PatchFrame (const PatchFrame& rhs)
: ComponentFrame(rhs)
{
    itsPatch = rhs.itsPatch ? new GeometryPatch(*rhs.itsPatch) : 0;
    SetGeometry(itsPatch);
}

inline PatchFrame::~PatchFrame()
{
    if(itsPatch)
        delete itsPatch;
}

#endif
