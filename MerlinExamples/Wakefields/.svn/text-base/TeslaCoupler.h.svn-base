// class TeslaCoupler
//
// coupler wake fields: - talk and paper(sign and decimal points errors removed)
//                        I. Zagordnov and M. Dohlus
//                        ILC/LCWS 07, Hamburg
//                      - steady state solution:
//                        M.Dohlus, I.Zagorodnov, E.Gjonaj and T.Weiland, 
//                        EPAC08, Genua, MOPP013 
// + tesla wake fields
//
// functions are selected by define statements for better performance

#ifndef _H_TeslaCoupler
#define _H_TeslaCoupler

#include "AcceleratorModel/CombinedWakeRF.h"
#include "NumericalUtils/Complex.h"
#include "EuclideanGeometry/Transform2D.h"
#include "NumericalUtils/PhysicalConstants.h"

using namespace PhysicalConstants;
using namespace std;

// Wxy         - coupler wakefield
// RF          - RF kick
// 
// old         - TDR coupler orientation
// new         - supposed modification of relative coupler orientation, see ILC/LCWS 07 talk
// oldRescaled - rescaled by 0.107 to model steady state solution, see MOPP013
// dummy       - contribution switched off

// functions are selected by define statements for better performance
#define oldWxy Wxy
//#define oldRescaledWxy Wxy
//#define newWxy Wxy
//#define dummyWxy Wxy
//----
//#define oldCouplerRFKick CouplerRFKick
//#define newCouplerRFKick CouplerRFKick
#define dummyCouplerRFKick CouplerRFKick

class TeslaCoupler : public CombinedWakeRF {
public:

	TeslaCoupler() : fac(2e12) {};
	// factor to calculate wake field kick
	// 2 because of kt -> constant ( = purely capacitive) wake
	// kV/nC -> V/C=V/nC*1e9*1e3

	// cavity wake fields - tesla wake fields
	//
	double Wlong(double z)  const {
		return 38.1e+12*(1.165*exp(-sqrt(z/3.65e-3))-0.165);
	};
	double Wtrans(double z) const {

	//	return (1290.0e+12*sqrt(z)-2600*z);       
	//	return 1.0e+12*(1290.0*sqrt(z)-2600.0*z);  //  V/C/m/module (per module)
		double arZ=sqrt(z/0.92e-03);               //  TESLA Report 2003-19
		return 1.21e14*(1.0-(1.0+arZ)*exp(-arZ));  //  V/C/m/m (per active length)
							   //  Merlin WakefielsProcess
							   //  does not use component length!
	};

	// coupler wake fields
	//
	// we need x,y since this is not just a transverse (dipole) wake field
	//
	// old upstream + downstream coupler design (sum of up+down from IZ&MD talk)
//	virtual Vector2D Wxy(double x, double y) const { // kV/nC 
	Vector2D oldWxy(double x, double y) const { // kV/nC 
		return fac*Vector2D( -0.021 + x*4.3 + y*0.07, -0.019 + x*0.03 - y*0.9 );
	};
	// old upstream + downstream coupler 
	// rescaled at (0,0) to match EPAC08, Genua, MOPP013
	Vector2D oldRescaledWxy(double x, double y) const { // kV/nC 
		return fac*0.107*Vector2D( -0.021 + x*4.3 + y*0.07, -0.019 + x*0.03 - y*0.9 );
	};
	// new downstream coupler design (sum of up+down from IZ&MD talk)
	Vector2D newWxy(double x,double y) const {
		return fac*Vector2D( -0.0025 + x*2.33 + y*0.04, -0.0002 - x*0.02 + y*1.1 );
	};
	Vector2D dummyWxy(double x, double y) const { 
		return Vector2D(0,0);
	} 

	// coupler RF kicks
	//
	// scaled kick = Re{Vt/Vz*exp(i*phi)} for particle a t - Vz=V_cavity, phi=phi0+2*pi*f*(t-t0)
	// a phi > means later than t0 - opposite sign to Merlin TWRFStructure::GetPhase()!
	//
	// old design
	Vector2D oldCouplerRFKickDown(double x, double y,double phi)  const {
		Complex a = polar(1., phi);

		// kap = Vt/Vz for downstream coupler
		// typo 0.00037 instead of 0.000037 in paper
		Complex kap_x(-0.000025 - x*.004  + y*.0029, 0.000052 - x*0.002 + y*0.00037);
		Complex kap_y( 0.000032 + x*.0029 + y*.0038, 0.000005 + x*.0005 + y*0.0018);

		return Vector2D(real(kap_x*a),real(kap_y*a));
		// komplex arithmetic is a bit slow

	}
	Vector2D oldCouplerRFKickUp(double x, double y,double phi)  const {
		Complex a = polar(1., phi);

		// kap = Vt/Vz for upstream coupler
		Complex kap_x=Complex(-0.000057 + x*0.0011 + y*0.0034, 0.000007 - x*0.0007 + y*0.00015);
		// typo in paper: 0.41+0.03i instead of the correct -0.41-0.03i
		Complex kap_y=Complex(-0.000041 + x*0.0034 - y*0.001,  -0.000003 + x*0.0002 + y*0.0006);

		return Vector2D(real(kap_x*a),real(kap_y*a));
		// komplex arithmetic is a bit slow
	}
	// new design - upstream unchanged
	Vector2D newCouplerRFKickUp(double x, double y,double phi) const {
		return oldCouplerRFKickUp(x,y,phi);
	};
	// new design - no numbers given by IZ&MD
	// we approx by mirroring downstream fields 
	Vector2D newCouplerRFKickDown(double x, double y,double phi) const {
		Complex a = polar(1., phi);

		// kap = Vt/Vz for downstream coupler
		Complex kap_x(-0.000025 - x*.004  + y*.0029,  0.000052 - x*0.002 + y*0.00037);
		// mirror i.e. - sign in y
		Complex kap_y(-0.000032 - x*.0029 - y*.0038, -0.000005 - x*.0005 - y*0.0018);

		return Vector2D(real(kap_x*a),real(kap_y*a));
	};

	Vector2D oldCouplerRFKick(double x, double y,double phi)  const {
		Complex a = polar(1., phi);
		// up+down old
		Complex kap_x(-0.000082  - x*0.0029  + y*0.0063, 0.000058  - x*0.0027 + y*0.00051);
		Complex kap_y(-0.0000092 + x*0.0063  + y*0.0028, 0.0000018 + x*0.0007 + y*0.0024);
		return Vector2D(real(kap_x*a),real(kap_y*a));
		// komplex arithmetic is a bit slow
	}
	Vector2D newCouplerRFKick(double x, double y,double phi)  const {
		Complex a = polar(1., phi);
		// up+down new
		Complex kap_x(-0.000082  - x*0.0029  + y*0.0063,  0.000058  - x*0.0027  + y*0.00051);
		Complex kap_y(-0.000074  + x*0.00049 - y*0.0048, -0.0000087 - x*0.00029 - y*0.0012);
		return Vector2D(real(kap_x*a),real(kap_y*a));
		// komplex arithmetic is a bit slow
	};

        Vector2D dummyCouplerRFKick(double x, double y,double phi)  const {
		return Vector2D(0,0);
        }	

private:
	double fac;
};

#endif
