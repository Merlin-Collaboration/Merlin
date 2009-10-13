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

#ifndef CorrectorWinding_h
#define CorrectorWinding_h 1

#include "merlin_config.h"
// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"
// ModelElement
#include "AcceleratorModel/ModelElement.h"

//	Represents a dipole corrector winding that can be added
//	to any RectMultipole component. The Bx and By field
//	values are in integrated strengths (Tesla.meter).

class CorrectorWinding : public ModelElement
{
public:
    CorrectorWinding (RectMultipole& aMagnet);

    void SetBx (double value);
    void SetBy (double value);
    double GetBx () const;
    double GetBy () const;

    //	Return the name of the element.
    virtual const string& GetName () const;

    //	Return the type string for the element.
    virtual const string& GetType () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

	//  Get the uniques beamline index for this frame
	size_t GetBeamlineIndex() const;
	void AppendBeamlineIndecies(std::vector<size_t>&) const;

private:

    RectMultipole* magnet;
};

inline const string& CorrectorWinding::GetName () const
{
    return magnet->GetName();
}

inline size_t CorrectorWinding::GetBeamlineIndex() const
{
	return magnet->GetBeamlineIndex();
}

inline void CorrectorWinding::AppendBeamlineIndecies(std::vector<size_t>& ivec) const
{
	magnet->AppendBeamlineIndecies(ivec);
}

#endif
