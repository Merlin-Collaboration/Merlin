/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef SupportStructure_h
#define SupportStructure_h 1

#include "merlin_config.h"
// AcceleratorSupport
#include "AcceleratorModel/Supports/AcceleratorSupport.h"
// SequenceFrame
#include "AcceleratorModel/Frames/SequenceFrame.h"
// Transform3D
#include "EuclideanGeometry/Transform3D.h"
// Rotation3D
#include "EuclideanGeometry/Rotation3D.h"

//	A SupportStructure represents a mechanical support
//	system, on which accelerator components can be mounted.
//	A SupportStructure is mounted to the ground by either
//	one Support placed at the centre of enclosed geometry,
//	or by two supports at the exit and entrance points.
//
//	SupportStructure and its associated Support class are
//	primarilly intended for application of ground motion and
//	vibration.

class SupportStructure : public SequenceFrame
{
public:

    typedef enum {simple,girder} Type;

    //	Copy constructor.
    SupportStructure (const SupportStructure& rhs);

    //	Destructor.
    ~SupportStructure ();

    //	Appends this structures support(s) to the SupportList.
    //	Returns the number of supports appended.
    int ExportSupports (AcceleratorSupportList& supports);

    //	Returns the local frame transformation. The result
    //	includes results of effects of the support offsets, as
    //	well as local transformations of the girder.
    virtual Transform3D GetLocalFrameTransform () const;

    //	When called, SupportStructure sets up is Accelerator
    //	Structure objects. This function should only be called
    //	after the AcceleratorModel is complete.
    virtual void ConsolidateConstruction ();

protected:

    SupportStructure (const string& id, Type type);

private:

    AcceleratorSupport* sup1;
    AcceleratorSupport* sup2;

    //	Updates (if necessary) the local frame transformation
    //	due to the support offsets.
    void UpdateSupportTransform () const;

    //	Rotation used to convert the support motion into the
    //	local entrance plane reference frame.
    Rotation3D Rg;

    //	Cached transformation state. Used to calculate the local
    //	frame transformation. Is recalculated if an offset of an
    //	accelerator support has changed.
    mutable Transform3D Ts;
};

//	Represents a long mount structure (girder) which has two
//	supports at either end. A girder can tilt under the
//	action of the two supports.

class GirderMount : public SupportStructure
{
public:

    explicit GirderMount (const std::string& id);

    //	Returns "SupportStructure".
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;
};

//	An accelerator mount which has a single support located
//	at the centre of the geometry (local origin).

class SimpleMount : public SupportStructure
{
public:

    explicit SimpleMount (const std::string& id);

    //	Returns "SupportStructure".
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;
};

inline GirderMount::GirderMount (const std::string& id)
        : SupportStructure(id,girder)
{}

inline SimpleMount::SimpleMount (const std::string& id)
        : SupportStructure(id,simple)
{}

#endif
