/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:29:37 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef LatticeFrame_h
#define LatticeFrame_h 1

#include "merlin_config.h"
// ModelElement
#include "AcceleratorModel/ModelElement.h"
// AcceleratorGeometry
#include "AcceleratorModel/AcceleratorGeometry.h"
// Transformable
#include "EuclideanGeometry/Transformable.h"

#define GLOBAL_FRAME (LatticeFrame*)0

// Class FrameTraverser
// An abstract iterator class that is used to
// traverse in order the nested frame hierachy

class LatticeFrame;

class FrameTraverser {
public:
    virtual ~FrameTraverser() {}
    virtual void ActOn(LatticeFrame* frame) =0;
};

//	A LatticeFrame is a ModelElement that provides a
//	coordinate system for a subsection of the accelerator
//	lattice. A lattice frame is located on the accelerator
//	geometry at a specified position s. LatticeFrames can be
//	nested in a frame hierachy, but may not partially
//	overlap. A LatticeFrame may represent the frame for a
//	single component (see class ComponentFrame), or may be
//	constructed from a contiguous sequence of LatticeFrames
//	(sub-frames).
//
//	Since each LatticeFrame contains a section of the
//	lattice, it represents the AcceleratorGeometry for that
//	section (which forms a sub-geometry of the total lattice
//	geometry). A LatticeFrame can therefore provide all the
//	geometry transformation functions within the its local
//	s-frame. Sub-frames of a LatticeFrame are located at
//	specific (local) s positions on that geometry.
//
//	LatticeFrames can be rotated and translated.
//	Rotating/translating a LatticeFrame rotates/translates
//	all the associated sub-frames. LatticeFrame provides
//	methods for accessing the various coordinate
//	transformations between various frames in the hierachy.
//	Two types of coordinate transformations are
//	distinguished:
//
//	- Frame Transformations: if the location of a Lattice
//	Frame on a super-frames geometry is s, then the frame
//	transformation to that super-frame is defined as the
//	coordinate transformation between the LatticeFrame's
//	local coordinate frame, and the coordinate frame defined
//	on the super-frame's geometry at s. Note that if no
//	rotations have been applied to any of the frames in the
//	hierachy (up to the specified super-frame), then the
//	frame transformation is always the identity.
//
//	- Physical Transformations: The physical transformation
//	from a LatticeFrame to a specified super-frame is the
//	total coordinate transformation between the local frame
//	of the LatticeFrame, and the local frame of the
//	super-frame. This transformation represents the physical
//	3-D location and orientation of the LatticeFrame's local
//	frame in the local frame of the super-frame. In general
//	it is not the identity, even when no transformations
//	have been applied.
//
//	Frame Transformations are most important during typical
//	beam dynamics tracking operations, while Physical
//	Transformations can be used for surveying.

class LatticeFrame : public ModelElement, public Transformable
{
public:

    typedef AcceleratorGeometry::BoundaryPlane BoundaryPlane;

    //	Constructor
    explicit LatticeFrame (const string& id = "");

    //	Copy constructor.
    LatticeFrame (const LatticeFrame& rhs);

    virtual ~LatticeFrame ();

    //	Frame Transformations.
    //
    //	Returns the frame transformation from the specified
    //	super-frame to this frame.
    virtual Transform3D GetFrameTransform (const LatticeFrame* sframe = GLOBAL_FRAME) const;

    //	Returns the frame transformation from the immediate
    //	super-frame, i.e. the result of all transformations
    //	applied to this frame.
    virtual Transform3D GetLocalFrameTransform () const;

    //	Physical Transformations.
    //
    //	Returns the physical transformation from the specified
    //	super-frame to this frame.
    Transform3D GetPhysicalTransform (const LatticeFrame* sframe = GLOBAL_FRAME) const;

    //	Returns the physical transformation from the immediate
    //	super-frame.
    Transform3D GetLocalPhysicalTransform () const;

    //	AcceleratorGeometry related operations.
    //
    //	Return the associated AcceleratorGeometry.
    const AcceleratorGeometry* GetGeometry () const;

    //	Returns the position of this frame's origin in the
    //	specified super-frames s-frame.
    double GetPosition (const LatticeFrame* sframe = GLOBAL_FRAME) const;

    //	Returns the position of this frame's origin in the
    //	immediate super-frames s-frame.
    double GetLocalPosition () const;

    //	Returns the AcceleratorGeometry transformation from s0
    //	to s (local s-frame).
    Transform3D GetGeometryTransform (double s0, double s) const throw (AcceleratorGeometry::BeyondExtent);

    //	Returns the AcceleratorGeometry transformation from the
    //	local origin to s (local s-frame).
    Transform3D GetGeometryTransform (double s) const throw (AcceleratorGeometry::BeyondExtent);

    //	Returns the transformation from the origin to the
    //	specified boundary plane.
    Transform3D GetGeometryTransform (BoundaryPlane p) const;

    //	Returns the transformation from the entrance plane frame
    //	to the exit plane frame.
    Transform3D GetTotalGeometryTransform () const;

    //	Returns the extent of this frame's geometry in the
    //	s-frame of the specified super-frame.
    AcceleratorGeometry::Extent GetGeometryExtent (const LatticeFrame* sframe = GLOBAL_FRAME) const;

    //	Returns the local extent of this frame's geometry, i.e.
    //	the extent in the local s-frame.
    AcceleratorGeometry::Extent GetLocalGeometryExtent () const;

    //	Returns the total arc-length of the associated geometry.
    double GetGeometryLength () const;

    //	Returns transformation to the specified boundary plane
    //	reference frame. The transformation includes all effects
    //	of nested boundary planes which occur at the plane.
    virtual Transform3D GetBoundaryPlaneTransform (BoundaryPlane p) const;

    //	Used during tracking. Returns the required
    //	transformation to the local entrance plane coordinate
    //	frame. Includes the effects of all nested frame
    //	boundaries which occur at this location.
    Transform3D GetEntrancePlaneTransform () const;

    //	Used during tracking. Returns the required
    //	transformation from the local exit plane coordinate
    //	frame. Includes the effects of all nested frame
    //	boundaries which occur at this location.
    Transform3D GetExitPlaneTransform () const;
/*
    //	Returns true if this frame has been locally transformed.
    bool IsTransformed () const;

    //	Translates the frame by the relative vector (dx,dy,dz).
    void Translate (double dx, double dy, double dz);

    //	Translates the frame by the relative vector X.
    void Translate (const Vector3D& X);

    //	Translates the frame along the current x-axis by dx.
    void TranslateX (double dx);

    //	Translates the frame along the current y-axis by dy.
    void TranslateY (double dy);

    //	Translates the frame along the current z-axis by dz.
    void TranslateZ (double dz);

    //	Rotates the frame about the current x-axis by angle.
    void RotateX (double angle);

    //	Rotates the frame about the current y-axis by angle.
    void RotateY (double angle);

    //	Rotates the frame about the current z-axis by angle.
    void RotateZ (double angle);
*/
    //	Transform the frame (with respect to the current axes)
    //	by the transformation t.
    void ApplyLocalFrameTransform (const Transform3D& t1);

    //	Set the local  frame transformation for this object.
    void SetLocalFrameTransform (const Transform3D& t);

    //	Clear the local frame transformation.
    void ClearLocalFrameTransform ();

    //	Function called after construction of the Accelerator
    //	Model is complete. Allows the nested frame hierachy to
    //	perform certain state checks and updates, which are only
    //	possible once the entire model is complete.
    virtual void ConsolidateConstruction ();

    //	Sets the position of the LatticeFrame on the immediate
    //	super-frame's geometry.
    void SetLocalPosition (double s);

    //	Set the super frame of this LatticeFrame object. Returns
    //	the old super frame.
    LatticeFrame* SetSuperFrame (LatticeFrame* aFrame);

    //	Replace subFrame with newSubFrame. Returns true if
    //	successful (i.e. subFrame is a sub-frame of this Lattice
    //	Frame).
    virtual bool ReplaceSubFrame (LatticeFrame* subFrame, LatticeFrame* newSubFrame);

    //	Returns true if this frame is the top-level frame of the
    //	hierachy.
    bool IsGlobalFrame () const;

    //	Returns the top-level (global)  frame of the hieracy.
    LatticeFrame* GetGlobalFrame () const;

    // Iterate a FrameTraverese object over the underlying
    // frame hierachy.
    virtual void Traverse(FrameTraverser& ft);

protected:

    //	Protected function called by concrete LatticeFrame
    //	classes during construction.
    void SetGeometry (const AcceleratorGeometry* geom);

    LatticeFrame* superFrame;

private:

    //	The s-position of this frame's origin on the super-frame
    //	s-frame.
    double s_0;

    //	The geometry associated with this frame.
    const AcceleratorGeometry *itsGeometry;

private:

    //	Returns true if aSubFrame is the first (entrance)
    //	sub-frame of this frame.
    virtual bool IsBoundaryPlane (BoundaryPlane p, const LatticeFrame* aSubFrame) const = 0;
    Transform3D LocalBoundaryPlaneTransform (BoundaryPlane p) const;
};

// Class LatticeFrame

inline LatticeFrame::LatticeFrame (const string& id)
        : ModelElement(id),Transformable(), s_0(0),superFrame(0),itsGeometry(0)
{}

inline LatticeFrame::LatticeFrame (const LatticeFrame& rhs)
        : ModelElement(rhs),Transformable(rhs),s_0(0),superFrame(0),itsGeometry(0)
{}

inline LatticeFrame::~LatticeFrame ()
{}

inline Transform3D LatticeFrame::GetLocalFrameTransform () const
{
    return local_T!=0 ? *local_T : Transform3D();
}

inline Transform3D LatticeFrame::GetLocalPhysicalTransform () const
{
    return GetPhysicalTransform(superFrame);
}

inline const AcceleratorGeometry* LatticeFrame::GetGeometry () const
{
    return itsGeometry;
}

inline double LatticeFrame::GetLocalPosition () const
{
    return s_0;
}

inline Transform3D LatticeFrame::GetGeometryTransform (double s0, double s) const throw (AcceleratorGeometry::BeyondExtent)
{
    if(itsGeometry!=0)
        return  itsGeometry->GetGeometryTransform(s0,s);
    else if(s0==s)
        return Transform3D();
    else
        throw AcceleratorGeometry::BeyondExtent();
}

inline Transform3D LatticeFrame::GetGeometryTransform (double s) const throw (AcceleratorGeometry::BeyondExtent)
{
    if(itsGeometry!=0)
        return  itsGeometry->GetGeometryTransform(s);
    else if(s==0)
        return Transform3D();
    else
        throw AcceleratorGeometry::BeyondExtent();
}

inline Transform3D LatticeFrame::GetGeometryTransform (BoundaryPlane p) const
{
    return (itsGeometry!=0)? itsGeometry->GetGeometryTransform(p) : Transform3D();
}

inline Transform3D LatticeFrame::GetTotalGeometryTransform () const
{
    return (itsGeometry!=0)? itsGeometry->GetTotalGeometryTransform() : Transform3D();
}

inline AcceleratorGeometry::Extent LatticeFrame::GetLocalGeometryExtent () const
{
    return (itsGeometry!=0)? itsGeometry->GetGeometryExtent() : AcceleratorGeometry::Extent(0,0);
}

inline double LatticeFrame::GetGeometryLength () const
{
    return (itsGeometry!=0) ? itsGeometry->GetGeometryLength() : 0;
}

inline Transform3D LatticeFrame::GetEntrancePlaneTransform () const
{
    return GetBoundaryPlaneTransform(AcceleratorGeometry::entrance);
}

inline Transform3D LatticeFrame::GetExitPlaneTransform () const
{
    return GetBoundaryPlaneTransform(AcceleratorGeometry::exit).inv();
}
/***
inline bool LatticeFrame::IsTransformed () const
{
    return local_T!=0 && !local_T->isIdentity();
}

inline void LatticeFrame::Translate (const Vector3D& X)
{
    Translate(X.x,X.y,X.z);
}

inline void LatticeFrame::TranslateX (double dx)
{
    Translate(dx,0,0);
}

inline void LatticeFrame::TranslateY (double dy)
{
    Translate(0,dy,0);
}

inline void LatticeFrame::TranslateZ (double dz)
{
    Translate(0,0,dz);
}

inline void LatticeFrame::Invalidate () const
{}
****/
inline void LatticeFrame::ConsolidateConstruction ()
{}

inline void LatticeFrame::SetLocalPosition (double s)
{
    s_0=s;
}

inline LatticeFrame* LatticeFrame::SetSuperFrame (LatticeFrame* aFrame)
{
    LatticeFrame* tmp=superFrame;
    superFrame = aFrame;
    return tmp;
}

inline bool LatticeFrame::ReplaceSubFrame (LatticeFrame* subFrame, LatticeFrame* newSubFrame)
{
    return false;
}

inline bool LatticeFrame::IsGlobalFrame () const
{
    return superFrame==0;
}

inline void LatticeFrame::SetGeometry (const AcceleratorGeometry* geom)
{
    itsGeometry=geom;
}

#endif
