/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _H_CombinedWakeRF
#define _H_CombinedWakeRF

#include <complex>
typedef std::complex<double> Complex;

#include "WakePotentials.h"
#include "Transform2D.h"
#include "PhysicalConstants.h"

using namespace PhysicalConstants;
using namespace std;

// Wxy            - Sum (up+downstream) of coupler wakefield
// CouplerRFKick  - Sum        "        of coupler RF kicks

/**
 * Interface class for e.g. coupler wakefields
 * i.e wakefields and RF kicks that depends on x,y
 * where wakefields depend on bunch charge
 * and RF kicks do not
 *
 * @see "../MerlinExamples/Wakefields" example 3
 * @see "../Merlin/CouplerWakeFieldProcess.cpp"
 */
class CombinedWakeRF: public WakePotentials
{
public:

	CombinedWakeRF()
	{
	}

	/**
	 * Coupler wake fields -- we need x,y since this is not just a transverse
	 * (dipole) wake field sum of upstream + downstream coupler
	 */
	virtual Vector2D Wxy(double x, double y) const = 0; // kV/nC

	/**
	 * Coupler RF kicks scaled kick equals
	 * \f[
	 *    \Re\left( V_t/V_z \times e^{i\phi} \right)
	 * \f]
	 * for particle at
	 * \f[
	 *    V_z=V_{cavity}, \quad \phi=\phi_0 + 2\pi f\times (t-t_0)
	 * \f]
	 * a larger \f$\phi\f$ means later than \f$t_0\f$ -- opposite sign to Merlin
	 * TWRFStructure::GetPhase
	 */
	virtual Vector2D CouplerRFKick(double x, double y, double phi) const = 0;
};

#endif
