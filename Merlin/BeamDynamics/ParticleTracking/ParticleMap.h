//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (1999)
* 
* file Merlin\BeamDynamics\ParticleTracking\ParticleMap.h
* last modified 08/10/01 02:14:16 PM
*/

/*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
* Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED. 
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


#ifndef ParticleMap_h
#define ParticleMap_h 1

#include "merlin_config.h"


namespace ParticleTracking {

class ParticleBunch;






//	An arbitrary map that can be applied to a ParticleBunch.







class ParticleMap
{
public:


    virtual ~ParticleMap ();




    //	Applies the map to the specified ParticleBunch.
    virtual ParticleBunch& Apply (ParticleBunch& bunch) const = 0;


    virtual void Invert () = 0;

protected:
private:
private:
};

// Class ParticleMap


inline ParticleMap::~ParticleMap ()
{


}


// Class ParticleMap

}; //end namespace ParticleTracking



#endif
