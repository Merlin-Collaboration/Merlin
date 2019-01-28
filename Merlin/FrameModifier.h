/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef FrameModifier_h
#define FrameModifier_h 1

#include "merlin_config.h"
#include "LatticeFrame.h"

/**
 *	A utility class which adds an additional layer to the
 *	frame hierarchy, and thus another level of coordinate
 *	transformations. FrameModifier effectively "wraps" a
 *	single LatticeFrame object. It can insert and remove
 *	itself in the lattice frame hierarchy.
 */
class FrameModifier: public LatticeFrame
{
public:

	/**
	 *	Constructor taking the LatticeFrame object above which
	 *	this is to be inserted.
	 */
	explicit FrameModifier(LatticeFrame* frame, const std::string& label = "");
	~FrameModifier();

	/**
	 *	Returns a copy of the sub-frame.
	 *	@return Copy of the sub-frame
	 */
	ModelElement* Copy() const;

	/**
	 *	Returns the type string of the sub-frame.
	 *	@return String type of the sub-frame
	 */
	const string& GetType() const;

	/**
	 *	Causes any cached state to be invalidated. The cached
	 *	state should be re-calculated if and when required.
	 */
	virtual void Invalidate() const;

	/**
	 *	Function called after construction of the Accelerator
	 *	Model is complete. Allows the nested frame hierarchy to
	 *	perform certain state checks and updates, which are only
	 *	possible once the entire model is complete.
	 */
	virtual void ConsolidateConstruction();

	/**
	 * Remove the modifier frame from the hierarchy
	 */
	void Remove();

private:

	LatticeFrame* subFrame;

	virtual bool IsBoundaryPlane(BoundaryPlane p, const LatticeFrame* aSubFrame) const;
};

#endif
