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

#ifndef PSvectorTransform3D_h
#define PSvectorTransform3D_h 1

#include "merlin_config.h"
// Transform3D
#include "EuclideanGeometry/Transform3D.h"
// PSTypes
#include "BeamModel/PSTypes.h"

//	Utility class for performing an arbitrary 3D coordinate
//	transformation on PSvector objects. The transformation
//	assumes that the particle(s) are fully relativistic.
//	Since small angle approximations are assumed, large
//	rotations about the x and y axis may lead to significant
//	errors.

class PSvectorTransform3D
{
public:

    PSvectorTransform3D (const Transform3D& tfrm);

    PSvector& Apply (PSvector& p) const;
    PSvectorArray& Apply (PSvectorArray& pv) const;
    PSvector& operator () (PSvector& p) const;

private:

    Transform3D T;
    bool bNoRot;
};

inline PSvector& PSvectorTransform3D::operator () (PSvector& p) const
{
    return Apply(p);
}

#endif
