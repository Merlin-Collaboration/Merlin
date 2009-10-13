// class CombinedWakePotentials
//
// interface class for e.g. coupler wakefields
// i.e wakefields and RF kicks that depends on x,y
// where wakefields depend on bunch charge
// and RF kicks do not; see MerlineExamples/Wakefields, example 3
// and BeamDynamics/ParticleTracking/CouplerWakeFieldProcess.h|cpp
// DK 3.12.08

#ifndef _H_CombinedWakeRF
#define _H_CombinedWakeRF

#include "AcceleratorModel/WakePotentials.h"
#include "NumericalUtils/Complex.h"
#include "EuclideanGeometry/Transform2D.h"
#include "NumericalUtils/PhysicalConstants.h"

using namespace PhysicalConstants;
using namespace std;

// Wxy            - Sum (up+downstream) of coupler wakefield
// CouplerRFKick  - Sum        "        of coupler RF kicks

class CombinedWakeRF : public WakePotentials {
public:

	CombinedWakeRF() {};
	// coupler wake fields
	//
	// we need x,y since this is not just a transverse (dipole) wake field
	//
	// sum of upstream + downstream coupler 
	virtual Vector2D Wxy(double x, double y) const = 0; // kV/nC 

	// coupler RF kicks
	//
	// scaled kick = Re{Vt/Vz*exp(i*phi)} for particle a t - Vz=V_cavity, phi=phi0+2*pi*f*(t-t0)
	// a phi > means later than t0 - opposite sign to Merlin TWRFStructure::GetPhase()!
	//
	virtual Vector2D CouplerRFKick(double x, double y,double phi)  const = 0;
};

#endif
