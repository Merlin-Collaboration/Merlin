/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef PatchFrame_h
#define PatchFrame_h 1

#include "merlin_config.h"
#include "ComponentFrame.h"
#include "GeometryPatch.h"

/**
 * A special ComponentFrame which represents a pure geometry patch.
 * There is no AcceleratorComponent associated with this frame.
 *
 * **NOTE** This his is a bit of a kludge, and there is probably a better
 *         way to do this!
 */

class PatchFrame: public ComponentFrame
{
public:

	/**
	 *	Constructor
	 */
	explicit PatchFrame(GeometryPatch* pg, const string& id = "");

	/**
	 *	Copy constructor.
	 */
	PatchFrame(const PatchFrame& rhs);

	/**
	 * Destructor
	 */
	virtual ~PatchFrame();
	virtual const string& GetType() const
	{
		_TYPESTR(PatchFrame)
	}

	virtual const string& GetName() const
	{
		return LatticeFrame::GetName();
	}

	virtual ModelElement* Copy() const
	{
		return new PatchFrame(*this);
	}

	virtual const Transform3D* GetEntranceGeometryPatch() const
	{
		return itsPatch ? itsPatch->GetTransformation() : nullptr;
	}

	void SetGeometryPatch(GeometryPatch* gp)
	{
		if(itsPatch)
		{
			delete itsPatch;
		}
		itsPatch = gp;
		SetGeometry(itsPatch);
	}

private:

	GeometryPatch* itsPatch;
};

inline PatchFrame::PatchFrame(GeometryPatch* pg, const string& id) :
	ComponentFrame(nullptr, id), itsPatch(pg)
{
	SetGeometry(itsPatch);
}

inline PatchFrame::PatchFrame(const PatchFrame& rhs) :
	ComponentFrame(rhs)
{
	itsPatch = rhs.itsPatch ? new GeometryPatch(*rhs.itsPatch) : nullptr;
	SetGeometry(itsPatch);
}

inline PatchFrame::~PatchFrame()
{
	if(itsPatch)
	{
		delete itsPatch;
	}
}

#endif
