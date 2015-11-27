#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
#include "AcceleratorModel/StdField/TWRFfield.h"
#include "AcceleratorModel/StdField/SWRFfield.h"
#include "BasicTransport/BasicTransportMaps.h"
#include "BasicTransport/MatrixMaps.h"
#include "BasicTransport/TransportMatrix.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"

#include "SymplecticIntegrators.h"
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
#include "BeamDynamics/ParticleTracking/Integrators/StdIntegrators.h"

#include <iostream>

namespace ParticleTracking{



namespace SYMPLECTIC {

DEF_INTG_SET(ParticleComponentTracker,StdISet)
ADD_INTG(DriftCI)
ADD_INTG(TWRFStructureCI)
ADD_INTG(SWRFStructureCI)
//~ ADD_INTG(THIN_LENS::SWRFStructureCI)
ADD_INTG(RectMultipoleCI)
ADD_INTG(SectorBendCI)
// from StdIntegrators
ADD_INTG(ParticleTracking::MarkerCI)
ADD_INTG(ParticleTracking::MonitorCI)
ADD_INTG(ParticleTracking::SolenoidCI)
END_INTG_SET

//template<> MAKE_DEF_INTG_SET(ParticleTracking::ParticleComponentTracker,ParticleTracking::SYMPLECTIC::StdISet)

	using namespace ParticleTracking;
	using namespace PhysicalUnits;
	using namespace PhysicalConstants;


	// Drift Map
	struct DriftMap {
	private:
		double ds;
	public:
		DriftMap(double _ds) : ds(_ds) {}
		void operator()(PSvector& v) const {
			double x0  = v.x();
			double y0  = v.y();
			double ct0 = v.ct();

			double d1  = 1.0 + v.dp();
			double k   = sqrt(d1*d1 - v.xp()*v.xp() - v.yp()*v.yp());

			v.x()  += v.xp()*ds/k;
			v.y()  += v.yp()*ds/k;
			v.ct() += ds - d1*ds/k;
		}
	};

	//Pole Face Rotation Map
	struct PoleFaceRotation {

	private:
		double R10;
		double R32;

	public:
		PoleFaceRotation(double h, SectorBend::PoleFace pf){
			const double theta = pf.rot;
			const double fint  = pf.fint;
			const double hgap  = pf.hgap;
			const double sinTheta = sin(theta);
			const double phi = 2.0*fint*hgap*h*(1+sinTheta*sinTheta)/cos(theta);
			R10 = h*tan(theta);
			R32 =-h*tan(theta-phi);
		}

		void operator()(PSvector& v) const {
			v.xp() += R10 * v.x();
			v.yp() += R32 * v.y();
		};
	};

	// Sector Bend Map (no quadrupole gradient)
	struct SectorBendMap {

	private:
		double h, ds;
	
	public:
		SectorBendMap(double _h, double _ds) : h(_h), ds(_ds) {}

		void operator()(PSvector& v) const {

			double& x0  = v.x();
			double& px0 = v.xp();
			double& y0  = v.y();
			double& py0 = v.yp();
			double& ct0 = v.ct();

			double  dp  = v.dp();
			double  d1  = 1.0 + dp;

			double  wx  = sqrt(h*h/d1);
			
			double  xs  = sin(  wx*ds);
			double  xc  = cos(  wx*ds);
			double  xs2 = sin(2*wx*ds);
			double  xc2 = cos(2*wx*ds);

			double  x1  =     x0*xc    + px0*xs*wx/h/h + dp*(1.0 - xc)/h;
			double  px1 =-h*h*x0*xs/wx + px0*xc        + dp*h*xs/wx;

			double  y1  = y0 + py0*ds/d1;
			double  py1 = py0;

			double  j2  = h*x0 - dp;

			double  c0  = -dp;
			double  c1  = -j2;
			double  c2  = -px0/h;
			double  c3  = -px0*px0/d1/d1/2.0;
			double  c4  =  px0*j2/h/d1;
			double  c5  = -j2*j2/d1/2.0;
			double  c6  = -py0*py0/d1/d1/2.0;

			double ct1  = ct0 + 
							(2*c0 + c3 + c5 + 2*c6)*ds/2.0 +
							c1*xs/wx + (c3 - c5)*xs2/wx/4.0 +
							c2*(1.0 - xc) + c4*(1.0 - xc2)/4.0;

			x0  = x1;
			px0 = px1;
			y0  = y1;
			py0 = py1;
			ct0 = ct1;
		}
	};

	// Sector Bend Map by Etienne Forest (no quadrupole gradient)
	struct SectorBendMapEF {

	private:
		double h, ds;
	
	public:
		SectorBendMapEF(double _h, double _ds) : h(_h), ds(_ds)	{}

		void operator()(PSvector& v) const {

			double& x0  = v.x();
			double& px0 = v.xp();
			double& y0  = v.y();
			double& py0 = v.yp();
			double& ct0 = v.ct();

			double  dp  = v.dp();
			double  d1  = 1.0 + dp;

			double  r  = 1.0/h;

			double  w  = sqrt(d1*d1-px0*px0-py0*py0) - h*(r + x0);
			double  rdpxds = -px0*sin(ds/r) + w*cos(ds/r);

			double  px1 = px0*cos(ds/r) + w*sin(ds/r);
			double  x1  = sqrt(d1*d1-px1*px1-py0*py0)/h - rdpxds/h - r;

			double  u  = asin(px0/sqrt(d1*d1-py0*py0)) - asin(px1/sqrt(d1*d1-py0*py0));

			double  y1  = y0 + py0*ds/h/r + py0*u/h;
			double  py1 = py0;

			double ct1  = ct0 + ds - d1*ds/h/r - d1*u/h;

			x0  = x1;
			px0 = px1;
			y0  = y1;
			py0 = py1;
			ct0 = ct1;
		}
	};

	// Sector Bend Map (with quadrupole gradient)
	struct CombinedFunctionSectorBendMap {

	private:
		double h, k1, ds;
	
	public:
		CombinedFunctionSectorBendMap(double _h, double _k1, double _ds) : h(_h), k1(_k1), ds(_ds) {}

		void operator()(PSvector& v) const {

			double xs,  xc,  ys,  yc;
			double xs2, xc2, ys2, yc2;

			double& x0  = v.x();
			double& px0 = v.xp();
			double& y0  = v.y();
			double& py0 = v.yp();
			double& ct0 = v.ct();

			double  dp  = v.dp();
			double  dp1 = 1.0 + v.dp();

			double  wx  = sqrt(fabs(h*h + k1)/dp1);
			double  wy  = sqrt(fabs(k1)/dp1);

			if((h*h + k1)>0) {
				xs  = sin(  wx*ds);
				xc  = cos(  wx*ds);
				xs2 = sin(2*wx*ds);
				xc2 = cos(2*wx*ds);
			} else {
				xs  = sinh(  wx*ds);
				xc  = cosh(  wx*ds);
				xs2 = sinh(2*wx*ds);
				xc2 = cosh(2*wx*ds);
			}

			if(k1>0) {
				ys  = sinh(  wy*ds);
				yc  = cosh(  wy*ds);
				ys2 = sinh(2*wy*ds);
				yc2 = cosh(2*wy*ds);
			} else {
				ys  = sin(  wy*ds);
				yc  = cos(  wy*ds);
				ys2 = sin(2*wy*ds);
				yc2 = cos(2*wy*ds);
			}

			double x1, px1;

			if((k1+h*h)==0) {
				x1  =          x0    +    px0*ds/(1.0+dp)        + dp*h*ds*ds/(1.0+dp)/2.0;
				px1 =                     px0*xc                 + dp*h*ds;
			} else {
				x1  =          x0*xc    + px0*xs*wx/fabs(k1+h*h) + dp*h*(1.0-xc)/(k1+h*h);
				px1 =-(k1+h*h)*x0*xs/wx + px0*xc                 + dp*h*xs/wx;
			}

			double y1  =     y0*yc    + py0*ys*wy/fabs(k1);
			double py1 =  k1*y0*ys/wy + py0*yc;

			double j1  = k1 + h*h;
			double j2  = (j1*x0 - h*dp);

			double c0  = -h*h*dp/j1;
			double c1  =  h*h*dp/j1 - h*x0;
			double c2  = -h*px0/j1;
			double c3  = -px0*px0/dp1/dp1/2.0;
			double c4  =  px0*j2/j1/dp1;
			double c5  = -j2*j2/j1/dp1/2.0;
			double c6  = -py0*py0/dp1/dp1/2.0;
			double c7  = -y0*py0/dp1;
			double c8  = -y0*y0*k1/dp1/2.0;
			
			double ct1 = ct0 + 
							(2*c0 + c3 + c5 + c6 - c8)*ds/2.0 +
							c1*xs/wx + (c3 - c5)*xs2/wx/4.0 + 
							c2*(1.0 - xc) + c4*(1.0 - xc2)/4.0 +
							(c6 + c8)*ys2/wy/4.0 - c7*(1.0 - yc2)/4.0;

			x0  = x1;
			px0 = px1;
			y0  = y1;
			py0 = py1;
			ct0 = ct1;
		}
	};

	// Quadrupole Map
	struct QuadrupoleMap {

	private:
		double k1, ds;

	public:
		QuadrupoleMap(double _k1, double _ds) : k1(_k1), ds(_ds) {}

		void operator()(PSvector& v) const {

			double& x0  = v.x();
			double& px0 = v.xp();
			double& y0  = v.y();
			double& py0 = v.yp();
			double& ct0 = v.ct();

			double dp1  = 1.0 + v.dp();
			double w    = sqrt(fabs(k1)/dp1);

			double xs,  xc,  ys,  yc;
			double xs2, xc2, ys2, yc2;

			if(k1>=0) {
				xs  = sin(w*ds);
				xc  = cos(w*ds);
				ys  = sinh(w*ds);
				yc  = cosh(w*ds);
				xs2 = sin(2*w*ds);
				xc2 = cos(2*w*ds);
				ys2 = sinh(2*w*ds);
				yc2 = cosh(2*w*ds);
			} else {
				xs  = sinh(w*ds);
				xc  = cosh(w*ds);
				ys  = sin(w*ds);
				yc  = cos(w*ds);
				xs2 = sinh(2*w*ds);
				xc2 = cosh(2*w*ds);
				ys2 = sin(2*w*ds);
				yc2 = cos(2*w*ds);
			}

			double x1  =     x0*xc   + px0*xs*w/fabs(k1);
			double px1 = -k1*x0*xs/w + px0*xc;

			double y1  =     y0*yc   + py0*ys*w/fabs(k1);
			double py1 =  k1*y0*ys/w + py0*yc;

			double c3  = -px0*px0/dp1/dp1/2.0;
			double c4  =  x0*px0/dp1;
			double c5  = -x0*x0*k1/dp1/2.0;
			double c6  = -py0*py0/dp1/dp1/2.0;
			double c7  = -y0*py0/dp1;
			double c8  = -y0*y0*k1/dp1/2.0;

			double ct1 = ct0 +
							(c3 + c5 + c6 - c8)*ds/2.0 +
							(c3 - c5)*xs2/w/4.0 + 
							c4*(1.0 - xc2)/4.0 +
							(c6 + c8)*ys2/w/4.0 - c7*(1.0 - yc2)/4.0;

			x0  = x1;
			px0 = px1;
			y0  = y1;
			py0 = py1;
			ct0 = ct1;
		}
	};

	// Thin Rectangular Multipole Map
	struct MultipoleKick {
	private:
		const MultipoleField& field;
		Complex scale;
		
	public:
		MultipoleKick(const MultipoleField& f, double ds, double P0, double phi=0) : field(f) {
			scale = ds*eV*SpeedOfLight/P0*Complex(cos(phi),sin(phi));
		}
		
		void operator()(PSvector& v) {
			double x=v.x();
			double y=v.y();
			Complex F = scale*field.GetField2D(x,y);
			v.xp() += -F.real();
			v.yp() +=  F.imag();
		}
	};

	//RF Structure Map
	struct RFStructureMap {
	private:
		double Vn, Ve, k, phi, dphi, cosPhi, d0;
		RMtrx& RM;
		bool fullacc;
	public:
		RFStructureMap(double Vnorm, double Verr, double kval, double phase, double phaseErr, RMtrx& R, bool full_accel)
			: Vn(Vnorm),Ve(Verr),k(kval),phi(phase),dphi(phaseErr),RM(R),fullacc(full_accel) {

			cosPhi = cos(phi);
			d0     = 1 + Vn*cosPhi;
		};
		void operator()(PSvector& v) const {
			RM.Apply(v);
			if(fullacc) {
				double k0  = sqrt((1.0+v.dp())*(1.0+v.dp()) - v.xp()*v.xp() - v.yp()*v.yp());
				v.xp() /= k0;
				v.yp() /= k0;
				v.dp() = (v.dp()+Vn*(1.0+Ve)*cos(phi+dphi-k*v.ct())-Vn*cosPhi)/d0;
				k0  = sqrt((1.0+v.dp())*(1.0+v.dp()) - v.xp()*v.xp() - v.yp()*v.yp());
				v.xp() *= k0;
				v.yp() *= k0;
			} else {
				v.dp() += Vn*(1.0+Ve)*cos(phi+dphi-k*v.ct());
			}
		}
	};

	//SWRF Structure Map
	struct SWRFStructureMap {
	private:
		double Vn, Ve, k, phi, dphi, cosPhi, VncosPhi, d0, lnd0, len;
	public:
		SWRFStructureMap(double Vnorm, double Verr, double kval, double phase, double phaseErr, double length)
			: Vn(Vnorm),Ve(Verr),k(kval),phi(phase),dphi(phaseErr),len(length) {

			cosPhi   = cos(phi);
			VncosPhi = Vn*cosPhi;
			d0       = 1 + VncosPhi;
			lnd0     = log(d0);
		};
		void operator()(PSvector& v) const {

			double k0  = sqrt((1.0+v.dp())*(1.0+v.dp()) - v.xp()*v.xp() - v.yp()*v.yp());
			double x0  = v.x();
			double xp0 = v.xp()/k0;
			double y0  = v.y();
			double yp0 = v.yp()/k0;

			double x1   =  (1.0 - 0.5*lnd0)*x0         + len*lnd0*xp0/VncosPhi;
			double xp1  = -d0*VncosPhi*lnd0*x0/4.0/len + (1.0 + 0.5*lnd0)*xp0/d0;

			double y1   =  (1.0 - 0.5*lnd0)*y0         + len*lnd0*yp0/VncosPhi;
			double yp1  = -d0*VncosPhi*lnd0*y0/4.0/len + (1.0 + 0.5*lnd0)*yp0/d0;

			v.dp()		= (v.dp()+Vn*(1.0+Ve)*cos(phi+dphi-k*v.ct())-Vn*cosPhi)/d0;

			double k1   = sqrt((1.0+v.dp())*(1.0+v.dp()) - xp1*xp1 - yp1*yp1);

			v.x()  = x1;
			v.xp() = xp1*k1;
			v.y()  = y1;
			v.yp() = yp1*k1;

		}
	};

	//Rosenzweig-Serafini RF Structure Map
	struct RSRFStructureMap {
	private:
		double Vn, Ve, k, phi, dphi, len, VncosPhi;
		double m11, m12, m21, m22;
	public:
		RSRFStructureMap(double Vnorm, double Verr, double kval, double phase, double phaseErr, double length)
			: Vn(Vnorm),Ve(Verr),k(kval),phi(phase),dphi(phaseErr),len(length) {

			double rt2      = sqrt(2.0);
			double rt8      = 2.0*rt2;

			double cosphi   = cos(phi);
			double sinphi   = sin(phi);

			VncosPhi        = Vn*cosphi;
			double d0       = 1.0 + VncosPhi;

			double alpha    = log(d0)/cosphi/rt8;
			double cosalpha = cos(alpha);
			double sinalpha = sin(alpha);

			m11  = cosalpha - rt2*cosphi*sinalpha;
			m12  = rt8*len*cosphi*sinalpha/VncosPhi;
			m21  = -(1.0 - 1.0/d0)*(sinphi/rt2 + 1.0/cosphi/rt8)*sinalpha/len;
			m22  = (cosalpha + rt2*cosphi*sinalpha)/d0;

			//~ cout<<std::setw(14)<<len;
			//~ cout<<std::setw(14)<<d0-1;
			//~ cout<<std::setw(14)<<alpha;
			//~ cout<<std::setw(14)<<m11;
			//~ cout<<std::setw(14)<<m12;
			//~ cout<<std::setw(14)<<m21;
			//~ cout<<std::setw(14)<<m22<<endl;

		};
		void operator()(PSvector& v) const {

			double k0  = sqrt((1.0+v.dp())*(1.0+v.dp()) - v.xp()*v.xp() - v.yp()*v.yp());
			double x0  = v.x();
			double xp0 = v.xp()/k0;
			double y0  = v.y();
			double yp0 = v.yp()/k0;

			double x1  = m11*x0 + m12*xp0;
			double xp1 = m21*x0 + m22*xp0;

			double y1  = m11*y0 + m12*yp0;
			double yp1 = m21*y0 + m22*yp0;

			v.dp()	   = (v.dp()+Vn*(1.0+Ve)*cos(phi+dphi-k*v.ct())-VncosPhi)/(1.0 + VncosPhi);

			double k1  = sqrt((1.0+v.dp())*(1.0+v.dp()) - xp1*xp1 - yp1*yp1);

			v.x()  = x1;
			v.xp() = xp1*k1;
			v.y()  = y1;
			v.yp() = yp1*k1;

		}
	};

	//Simple RF Structure Map
	//Does not include correct transverse focusing
	struct SimpleRFStructureMap {
	private:
		double Vn, Ve, k, phi, dphi, cosPhi, VncosPhi, d0, len;
	public:
		SimpleRFStructureMap(double Vnorm, double Verr, double kval, double phase, double phaseErr, double length)
			: Vn(Vnorm),Ve(Verr),k(kval),phi(phase),dphi(phaseErr),len(length) {

			cosPhi   = cos(phi);
			VncosPhi = Vn*cosPhi;
			d0       = 1 + VncosPhi;
		};
		void operator()(PSvector& v) const {

			double k0  = sqrt((1.0+v.dp())*(1.0+v.dp()) - v.xp()*v.xp() - v.yp()*v.yp());
			double x0  = v.x();
			double xp0 = v.xp()/k0;
			double y0  = v.y();
			double yp0 = v.yp()/k0;

			double x1   = x0 + xp0*len; 
			double xp1  = xp0/(1.0+VncosPhi);

			double y1   = y0 + yp0*len;
			double yp1  = yp0/(1.0+VncosPhi);

			v.dp()		= (v.dp()+Vn*(cos(phi-k*v.ct())-cosPhi))/d0;

			double k1   = sqrt((1.0+v.dp())*(1.0+v.dp()) - xp1*xp1 - yp1*yp1);

			v.x()  = x1;
			v.xp() = xp1*k1;
			v.y()  = y1;
			v.yp() = yp1*k1;
		};
	};

	//RT Map
	struct ApplyRTMap {
	private:
		RTMap* m;

	public:
		ApplyRTMap(RTMap* M) : m(M) {};
		void operator()(PSvector& v) const {
			m->Apply(v);
		};
	};



	// Functors for applying maps to a bunch

	inline void ApplyDriftMap(ParticleBunch* bunch, double ds) {
		if(ds!=0)
			for_each(bunch->begin(),bunch->end(),DriftMap(ds));
	}

	inline void ApplyMultipoleKick(ParticleBunch* bunch, MultipoleField& field, double ds, double P0, double q) {
		if(ds!=0)
			for_each(bunch->begin(),bunch->end(),MultipoleKick(field, ds, P0));
	}

	inline void ApplyPoleFaceRotation(ParticleBunch* bunch, double h, const SectorBend::PoleFace& pf) {
		for_each(bunch->begin(),bunch->end(),PoleFaceRotation(h,pf));
	}

	inline void ApplySectorBendMap(ParticleBunch* bunch, double h, double ds) {
		if(ds!=0)
			if(h==0)
				for_each(bunch->begin(),bunch->end(),DriftMap(ds));
			else
				for_each(bunch->begin(),bunch->end(),SectorBendMap(h, ds));
	}

	inline void ApplyCombinedFunctionSectorBendMap(ParticleBunch* bunch, double h, double k1, double ds) {
		if(ds!=0)
			for_each(bunch->begin(),bunch->end(),CombinedFunctionSectorBendMap(h, k1, ds));
	}

	inline void ApplyQuadrupoleMap(ParticleBunch* bunch, double k1, double ds) {
		if(ds!=0)
			for_each(bunch->begin(),bunch->end(),QuadrupoleMap(k1, ds));
	}

	inline void ApplyRFStructureMap(ParticleBunch* bunch, double Vnorm, double Verr, double kval, double phase, double phaseErr, RMtrx& RM, bool full_accel) {
		for_each(bunch->begin(),bunch->end(),RFStructureMap(Vnorm, Verr, kval, phase, phaseErr, RM, full_accel));
	}

	inline void ApplySWRFStructureMap(ParticleBunch* bunch, double Vnorm, double Verr, double kval, double phase, double phaseErr, double length) {
		for_each(bunch->begin(),bunch->end(),RSRFStructureMap(Vnorm, Verr, kval, phase, phaseErr, length));
	}

	inline void ApplySimpleRFStructureMap(ParticleBunch* bunch, double Vnorm, double Verr, double kval, double phase, double phaseErr, double length) {
		for_each(bunch->begin(),bunch->end(),SimpleRFStructureMap(Vnorm, Verr, kval, phase, phaseErr, length));
	}



	
	
	// TrackStep Routines

	void DriftCI::TrackStep (double ds) {
		ApplyDriftMap(currentBunch, ds);
	}

	void SectorBendCI::TrackEntrance() {
		double h = (*currentComponent).GetGeometry().GetCurvature();
		const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
		if(pfi.entrance!=0)
			ApplyPoleFaceRotation(currentBunch, h, *pfi.entrance);
	}

	void SectorBendCI::TrackExit() {
		double h = (*currentComponent).GetGeometry().GetCurvature();
		const SectorBend::PoleFaceInfo& pfi = currentComponent->GetPoleFaceInfo();
		if(pfi.exit!=0)
			ApplyPoleFaceRotation(currentBunch, h, *pfi.exit);
	}

	void SectorBendCI::TrackStep (double ds) {

		double h = currentComponent->GetGeometry().GetCurvature();

		double tilt = currentComponent->GetGeometry().GetTilt();

		if(tilt!=0) {
			RMtrx Rr;
			TransportMatrix::Srot(tilt,Rr.R);
			Rr.Apply(currentBunch->GetParticles());
		}

		MultipoleField& field = currentComponent->GetField();
		const double P0   = currentBunch->GetReferenceMomentum();
		const double q    = currentBunch->GetChargeSign();
		const double Pref = currentComponent->GetMatchedMomentum(q);
		const double brho = P0/eV/SpeedOfLight;
		int np = field.HighestMultipole();

		assert(Pref>0);
		
		const Complex b0 = field.GetCoefficient(0);
		const Complex K1 = (np>0) ? q*field.GetKn(1,brho) : Complex(0);

		// We need to split the magnet for a kick if the following is true
		bool splitMagnet = b0.imag()!=0 || K1.imag()!=0 || np>1;
		double len = splitMagnet ? ds/2.0 : ds; 

		if(h==0 && K1.real()==0)
			ApplyDriftMap(currentBunch,len);

		if(h==0 && K1.real()!=0)
			ApplyQuadrupoleMap(currentBunch,K1.real(),len);

		if(h!=0 && K1.real()==0)
			ApplySectorBendMap(currentBunch,h,len);

		if(h!=0 && K1.real()!=0)
			ApplyCombinedFunctionSectorBendMap(currentBunch,h,K1.real(),len);


		// If we have split the magnet, we need to apply the kick approximation,
		// and then re-apply the map M
		if(splitMagnet) {

			// First we set the real parts of the dipole and quad fields to zero,
			// since these components have been modeled in the matrix
			Complex b1=field.GetCoefficient(1);
			field.SetCoefficient(0,Complex(0,b0.imag()));
			field.SetCoefficient(1,Complex(0,b1.imag()));

			// Apply the integrated kick, and then track through the linear second half
			ApplyMultipoleKick(currentBunch, field, ds, P0, q);

			if(h==0 && K1.real()==0)
				ApplyDriftMap(currentBunch,len);

			if(h==0 && K1.real()!=0)
				ApplyQuadrupoleMap(currentBunch,K1.real(),len);

			if(h!=0 && K1.real()==0)
				ApplySectorBendMap(currentBunch,h,len);

			if(h!=0 && K1.real()!=0)
				ApplyCombinedFunctionSectorBendMap(currentBunch,h,K1.real(),len);

			// Remember to set the components back
			field.SetCoefficient(0,b0);
			field.SetCoefficient(1,b1);
		}

		if(tilt!=0) {
			RMtrx Rr;
			TransportMatrix::Srot(-tilt,Rr.R);
			Rr.Apply(currentBunch->GetParticles());
		}
	}

	void RectMultipoleCI::TrackStep (double ds)
	{

		// Here we use a matrix to represent the quadrupole term, and a
		// single kick at the centre of the element for the other multipoles,
		// including any dipole term.
		using namespace TLAS;

		if(ds==0) return;

		double P0 = currentBunch->GetReferenceMomentum();
		double q = currentBunch->GetChargeSign();
		double brho = P0/eV/SpeedOfLight;

		MultipoleField& field = currentComponent->GetField();

		const Complex cK1 = q*field.GetKn(1,brho);
		bool splitMagnet = field.GetCoefficient(0)!=0.0 || field.HighestMultipole()>1;
		double len = splitMagnet ? ds/2 : ds;

		RdpMtrx M(2);

		if(cK1!=0.0) {
			double K1 = cK1.imag()==0 ? cK1.real() : abs(cK1);
			TransportMatrix::QuadrupoleR(len,K1,M.R);
			TransportMatrix::QuadrupoleT(len,K1,M.T);

			if(cK1.imag()!=0) {
				// Need to rotate the map
				double a = arg(cK1)/2;
				RealMatrix Rr(4,4);
				TransportMatrix::Srot(a,Rr);
				M.R = Rr*M.R*Transpose(Rr);
				M.T = Rr*M.T*Transpose(Rr);
			}

			M.Apply(currentBunch->GetParticles());
		}
		else
			ApplyDriftMap(currentBunch,len);

		if(splitMagnet) {
			Complex b1 = field.GetCoefficient(1);
			field.SetCoefficient(1,Complex(0));
			for_each(currentBunch->begin(),currentBunch->end(),MultipoleKick(field,ds,P0,q));
			if(b1!=0.0)
				M.Apply(currentBunch->GetParticles());
			else
				ApplyDriftMap(currentBunch,len);
			field.SetCoefficient(1,b1);
		}
	}	

	void TWRFStructureCI::TrackStep (double ds) {
		const TWRFfield& field = dynamic_cast<TWRFfield&>(currentComponent->GetField());
		double g = field.GetAmplitude();
		
		if(g==0) {
			ApplyDriftMap(currentBunch,ds);
			return;
		}

		double k    = field.GetK();
		double phi  = field.GetPhase();
		double dphi = 0;
		double dV   = 0;
		double P0   = currentBunch->GetReferenceMomentum();

		RMtrx RM(2);
		TransportMatrix::Drift(0,RM.R);

		ApplyDriftMap(currentBunch,ds/2.0);
		ApplyRFStructureMap(currentBunch,g*ds/P0,dV,k,phi,dphi,RM,false);
		ApplyDriftMap(currentBunch,ds/2.0);
	}

	void SWRFStructureCI::TrackStep (double ds) {
		const SWRFfield& field = dynamic_cast<SWRFfield&>(currentComponent->GetField());
		double g   = field.GetAmplitude();
		
		if(g==0) {
			ApplyDriftMap(currentBunch,ds);
			return;
		}

		double f   = field.GetFrequency();
		double phi = field.GetPhase();
		double dphi= 0;
		double dV  = 0;
		double P0  = currentBunch->GetReferenceMomentum();
		double k   = twoPi*f/SpeedOfLight;
		double lambdaOver2 = SpeedOfLight/f/2;

		int ncells = static_cast<int>(ds/lambdaOver2);
		assert(ncells*lambdaOver2 == ds);

		ApplySWRFStructureMap(currentBunch,g*ds/P0,dV,k,phi,dphi,ds);

		currentBunch->IncrReferenceMomentum(g*ds*cos(phi));
	}

} // end namespace SYMPLECTIC
}

