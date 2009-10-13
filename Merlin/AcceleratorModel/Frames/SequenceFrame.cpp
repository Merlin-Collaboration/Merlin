/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/28 19:36:24 $
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#include "NumericalUtils/utils.h"
// AcceleratorGeometry
#include "AcceleratorModel/AcceleratorGeometry.h"
// SequenceFrame
#include "AcceleratorModel/Frames/SequenceFrame.h"
// deleters
#include "stdext/deleters.h"

class SequenceGeometry;

//	Implementation of an AcceleratorGeometry class used by
//	SequenceFrame to construct the required geometry object
//	to pass to its LatticeFrame base class during
//	construction (via LatticeFrame::SetGeometry).

class SequenceGeometry : public AcceleratorGeometry
{
public:

    SequenceGeometry (SequenceFrame::FrameList* frames, SequenceFrame::Origin origin);

    //	Copy constructor.
    SequenceGeometry (const SequenceGeometry& rhs);

    ~SequenceGeometry ();

    //	Return the three-dimensional transformation between the
    //	frame at s0 and the frame at s. s and s0 are in the
    //	geometry's s-frame, and must be within the geometry
    //	extents.
    virtual Transform3D GetGeometryTransform (double s0, double s) const throw (BeyondExtent);

    //	Returns the transformation from the geometry origin to
    //	the specified boundary plane.
    virtual Transform3D GetGeometryTransform (BoundaryPlane p) const;

    //	Returns the transformation from the entrance plane frame
    //	to the exit plane frame.
    virtual Transform3D GetTotalGeometryTransform () const;

    //	Returns the local extent of this geometry.
    virtual AcceleratorGeometry::Extent GetGeometryExtent () const;

    //	Returns the total arc-length of the geometry.
    virtual double GetGeometryLength () const;

    //	Calculates the cached transform state.
    void CalculateCachedTransforms ();

    //	Calculates the cached transforms for centred origin.
    void CalculateCentreCT ();

    double len_t;
    SequenceFrame::Origin omode;

private:

    SequenceFrame::FrameList* theFrameList;

    struct TCache
    {
        Transform3D t_in;
        Transform3D t_out;
    };

    TCache* t_cache;
};

inline SequenceGeometry::SequenceGeometry (SequenceFrame::FrameList* frames, SequenceFrame::Origin origin)
        : len_t(0),omode(origin),t_cache(0),theFrameList(frames)
{}

SequenceFrame::SequenceFrame (const string& id, Origin originLoc)
        : LatticeFrame(id)
{
    itsSeqGeom = new SequenceGeometry(&subFrames,originLoc);
    SetGeometry(itsSeqGeom);
}

SequenceFrame::SequenceFrame (const SequenceFrame& rhs)
        : LatticeFrame(rhs)
{
    itsSeqGeom = new SequenceGeometry(*(rhs.itsSeqGeom));
    SetGeometry(itsSeqGeom);
    CopySubFrames(rhs.subFrames);
}

SequenceFrame::~SequenceFrame ()
{
    if(itsSeqGeom)
        delete itsSeqGeom;
}

void SequenceFrame::AppendFrame (LatticeFrame& aFrame)
{
    subFrames.push_back(&aFrame);
    aFrame.SetSuperFrame(this);

    // Update the geometry length
    itsSeqGeom->len_t += aFrame.GetGeometryLength();
}

void SequenceFrame::ConsolidateConstruction ()
{
    // First we need to set up the cached transform information
    // We assume that the origin of the sequence geometry is at
    // the arc centre.

	double s=0;
    switch(itsSeqGeom->omode) {
    case originAtEntrance:
        s=0;
        break;
    case originAtCenter:
        s=-(itsSeqGeom->len_t)/2;
        break;
    case originAtExit:
        s=-(itsSeqGeom->len_t);
        break;
    }

    for(FrameList::iterator fi = subFrames.begin(); fi!=subFrames.end(); ++fi) {
        // First update the local position of each sub-frame
        // to reflect the centred origin
        AcceleratorGeometry::Extent ext = (*fi)->GetLocalGeometryExtent();
        s-=ext.first;
        (*fi)->SetLocalPosition(s);
        s+=ext.second;

        // now forward consolidation to the sub-frame
        (*fi)->ConsolidateConstruction();
    }

    itsSeqGeom->CalculateCachedTransforms();
}

void SequenceFrame::Invalidate () const
{
    if(!subFrames.empty()){
        (subFrames.front())->Invalidate();
        if(subFrames.size()>1)
            (subFrames.back())->Invalidate();
    }
}

ModelElement* SequenceFrame::Copy () const
{
    return new SequenceFrame(*this);
}

bool SequenceFrame::ReplaceSubFrame (LatticeFrame* subFrame, LatticeFrame* newSubFrame)
{
    FrameList::iterator fi = std::find(subFrames.begin(),subFrames.end(),subFrame);
    if(fi!=subFrames.end()) {
        *fi = newSubFrame;
        newSubFrame->SetSuperFrame(this);
        return true;
    }
    else
        return false;
}

const string& SequenceFrame::GetType () const
{
    _TYPESTR(SequenceFrame);
}

void SequenceFrame::CopySubFrames (const list<LatticeFrame*>& frames)
{
    subFrames.clear();
    for(FrameList::const_iterator fi = frames.begin(); fi!=frames.end(); ++fi) {
        LatticeFrame* newCopy = static_cast<LatticeFrame*>((*fi)->Copy());
        subFrames.push_back(newCopy);
        newCopy->SetLocalPosition((*fi)->GetLocalPosition());
        newCopy->SetSuperFrame(this);
    }
}

bool SequenceFrame::IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const
{
    if(
        (p==AcceleratorGeometry::entrance && subFrames.front()==aSubFrame)||
        (p==AcceleratorGeometry::exit && subFrames.back()==aSubFrame))
        return true;
    else
        return false;
}

SequenceGeometry::SequenceGeometry (const SequenceGeometry& rhs)
{}

SequenceGeometry::~SequenceGeometry ()
{
    if(t_cache)
        delete t_cache;
}

Transform3D SequenceGeometry::GetGeometryTransform (double s0, double s) const throw (BeyondExtent)
{
    if(fequal(s,s0))
        return Transform3D();

    if(omode==SequenceFrame::originAtCenter) {
        s+=len_t/2;
        s0+=len_t/2;
    }
    else if(omode==SequenceFrame::originAtExit) {
        s+=len_t;
        s0+=len_t;
    }

    if(s<0||s>len_t||s0<0||s0>len_t)
        throw BeyondExtent();

	if(fequal(s0,0) && fequal(s,len_t))
		return GetTotalGeometryTransform();

    bool inv_t;
    if(inv_t=s0>s)
        std::swap(s,s0);

    SequenceFrame::FrameList::const_iterator fi = theFrameList->begin();
    double l=0;

    // First identify in which sub-frame s0 is
    for(;fi!=theFrameList->end();fi++) {
        l=(*fi)->GetGeometryLength();
        if(fequal(s0,l)||s0<l)
            break;
        s0-=l;
        s-=l;
    }

	assert(fi!=theFrameList->end());

    // s0 is now the distance from the entrance plane of *fi
    // to the requested starting point.
    // Now check if s is also within this frame. If so,
    // return the local transform
    Extent ext = (*fi)->GetLocalGeometryExtent();
    if(fequal(s,l)||s<l)
        return (*fi)->GetGeometryTransform(ext.first+s0,ext.first+s);

    Transform3D t = (*fi)->GetGeometryTransform(ext.first+s0,ext.second);
    s-=l;
    fi++;

    for(;fi!=theFrameList->end();fi++) {
        l=(*fi)->GetGeometryLength();
        if(fequal(s,l)||s<l)
            break;
        t=((*fi)->GetTotalGeometryTransform())*t;
        s-=l;
    }
    
	assert(fi!=theFrameList->end());

    // s is now the distance from the entrance plane of the *fi

    ext = (*fi)->GetLocalGeometryExtent();
    t = ((*fi)->GetGeometryTransform(ext.first,ext.first+s))*t;

    return inv_t ? t.inv() : t;
}

Transform3D SequenceGeometry::GetGeometryTransform (BoundaryPlane p) const
{
    return p==entrance ? t_cache->t_in : t_cache->t_out;
}

Transform3D SequenceGeometry::GetTotalGeometryTransform () const
{
    if(t_cache)
        return (t_cache->t_out)*((t_cache->t_in).inv());

    Transform3D t;
    for(SequenceFrame::FrameList::const_iterator fi = theFrameList->begin(); fi!=theFrameList->end(); fi++)
        t=((*fi)->GetTotalGeometryTransform())*t;

    return t;
}

AcceleratorGeometry::Extent SequenceGeometry::GetGeometryExtent () const
{
    Extent ext;
    switch(omode) {
    case SequenceFrame::originAtEntrance:
        ext.first=0;
        ext.second=len_t;
        break;
    case SequenceFrame::originAtCenter:
        ext.first=-len_t/2;
        ext.second=len_t/2;
        break;
    case SequenceFrame::originAtExit:
        ext.first=-len_t;
        ext.second=0;
        break;
    }

    return ext;
}

double SequenceGeometry::GetGeometryLength () const
{
    return len_t;
}

void SequenceGeometry::CalculateCachedTransforms ()
{
    // Set up cached transformations for efficiency
    if(t_cache==0)
        t_cache = new TCache;

    switch(omode) {
    case SequenceFrame::originAtEntrance:
        t_cache->t_out = GetTotalGeometryTransform();
        break;
    case SequenceFrame::originAtCenter:
        //		t_cache->t_in = GetGeometryTransform(0,-len_t/2);
        //		t_cache->t_out = GetGeometryTransform(0,len_t/2);
        CalculateCentreCT();
        break;
    case SequenceFrame::originAtExit:
        t_cache->t_in = GetTotalGeometryTransform().inv();
        break;
    }

}

void SequenceGeometry::CalculateCentreCT ()
{
    double s=0;
    double s0 = len_t/2;

    SequenceFrame::FrameList::const_iterator fi = theFrameList->begin();
    Transform3D t;

    for(;fi!=theFrameList->end();fi++) {
        s+=(*fi)->GetGeometryLength();
        if(fequal(s,s0)||s<s0)
            t=((*fi)->GetTotalGeometryTransform())*t;

        if(fequal(s,s0)) {
            t_cache->t_in = t.inv();
            t=Transform3D();
            break;
        }
        else if(s>s0) {
            Extent ext = (*fi)->GetLocalGeometryExtent();
            double s1 = ext.first+s-s0;
            t = (*fi)->GetGeometryTransform(ext.first,s1)*t;
            t_cache->t_in = t.inv();
            t = (*fi)->GetGeometryTransform(s1,ext.second);
            break;
        }
    }

    assert(fi!=theFrameList->end());
    fi++;

    for(;fi!=theFrameList->end();fi++)
        t=((*fi)->GetTotalGeometryTransform())*t;

    t_cache->t_out = t;
}


void SequenceFrame::Traverse(FrameTraverser &ft)
{
    // First act on this frame, and then
    // iterate the traverser over the sub-frames
    ft.ActOn(this);

    for(FrameList::iterator fi = subFrames.begin(); fi!=subFrames.end(); fi++)
        (*fi)->Traverse(ft);
}

void SequenceFrame::AppendBeamlineIndecies(std::vector<size_t>& ivec) const
{
    for(SequenceFrame::FrameList::const_iterator fi = subFrames.begin(); 
		fi!=subFrames.end(); fi++){
			(*fi)->AppendBeamlineIndecies(ivec);
	}
}


