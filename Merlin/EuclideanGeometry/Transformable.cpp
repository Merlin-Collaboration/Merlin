/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:17:37 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
// Transformable
#include "EuclideanGeometry/Transformable.h"

// macro for transformations
#define _TRNSFM(func) \
	if(local_T) (*local_T)*=Transform3D::func; \
	else local_T = new Transform3D(Transform3D::func); \
	Invalidate();


Transformable::~Transformable()
{
    if(local_T)
        delete local_T;
}

void Transformable::Translate (double dx, double dy, double dz)
{
    _TRNSFM(translation(dx,dy,dz));
}

void Transformable::RotateX (double angle)
{
    _TRNSFM(rotationX(angle));
}

void Transformable::RotateY (double angle)
{
    _TRNSFM(rotationY(angle));
}

void Transformable::RotateZ (double angle)
{
    _TRNSFM(rotationZ(angle));
}

void Transformable::ClearTransform ()
{
    if(local_T){
        delete local_T;
        local_T=0;
    }
    Invalidate();
}
