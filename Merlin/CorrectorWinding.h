/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CorrectorWinding_h
#define CorrectorWinding_h 1

#include "merlin_config.h"
#include "RectMultipole.h"
#include "ModelElement.h"

/**
 *	Represents a dipole corrector winding that can be added
 *	to any RectMultipole component. The Bx and By field
 *	values are in integrated strengths (Tesla.meter).
 */

class CorrectorWinding: public ModelElement
{
public:
	CorrectorWinding(RectMultipole& aMagnet);

	void SetBx(double value);
	void SetBy(double value);
	double GetBx() const;
	double GetBy() const;

	/**
	 *	Return the name of the element.
	 */
	virtual const string& GetName() const;

	/**
	 *	Return the type string for the element.
	 */
	virtual const string& GetType() const;

	/**
	 *	Virtual constructor.
	 */
	virtual ModelElement* Copy() const;

	/**
	 *  Get the unique beamline index for this frame
	 */
	size_t GetBeamlineIndex() const;
	void AppendBeamlineIndexes(std::vector<size_t>&) const;

private:

	RectMultipole* magnet;
};

inline const string& CorrectorWinding::GetName() const
{
	return magnet->GetName();
}

inline size_t CorrectorWinding::GetBeamlineIndex() const
{
	return magnet->GetBeamlineIndex();
}

inline void CorrectorWinding::AppendBeamlineIndexes(std::vector<size_t>& ivec) const
{
	magnet->AppendBeamlineIndexes(ivec);
}

#endif
