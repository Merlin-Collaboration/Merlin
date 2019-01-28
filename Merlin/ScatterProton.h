/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "PSvector.h"
#include "MADInterface.h"
#include "utils.h"
#include "TiltedAperture.hpp"

pair<double, double> CoulombScatterp(double x, double theta0);
int ScatterProton(PSvector& p, double x, double E0, const Aperture* tap);
//int ScatterProtonQ(PSvectorQ& p, double x, double E0,const  Aperture* tap);
