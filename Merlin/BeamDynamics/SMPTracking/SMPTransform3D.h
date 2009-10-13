/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (2000)
* 
* file Merlin\BasicTransport\PSvectorTransform3D.h
* last modified 04/04/01 14:41:34
*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
*
* Copyright (c) 2000 by The Merlin Collaboration.  
* ALL RIGHTS RESERVED. 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/

#ifndef SMPTransform3D_h
#define SMPTransform3D_h 1

#include "merlin_config.h"
#include "EuclideanGeometry/Transform3D.h"
#include "BeamDynamics/SMPTracking/SliceMacroParticle.h"
#include "BasicTransport/RMap.h"

namespace SMPTracking {

class SMPTransform3D {
public:

    SMPTransform3D(const Transform3D& tfrm);

    // Apply (approximate) transformation
    SliceMacroParticle& Apply (SliceMacroParticle& p) const;
    SliceMacroParticle& operator () (SliceMacroParticle& p) const{
        return Apply(p);
    }

private:

    R2Map R;
    double delta_x,delta_y,theta_x,theta_y;
    bool bNoRot;
    bool nullRotation;
};

}; // end namespace SMPTracking

#endif // SMPTransform3D_h
