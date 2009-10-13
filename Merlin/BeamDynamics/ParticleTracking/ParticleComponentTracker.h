/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\BeamDynamics\ParticleTracking\ParticleBunchTracker.h
 * last modified 09/11/01 04:10:33 PM
 */
/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 2001 by The Merlin Collaboration.
 * - ALL RIGHTS RESERVED - 
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

#ifndef ParticleBunchTracker_h
#define ParticleBunchTracker_h 1

#include "merlin_config.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

namespace ParticleTracking {
typedef TBunchCMPTracker<ParticleBunch> ParticleComponentTracker;
};


#endif
