#include "BeamDynamics/ParticleTracking/CouplerWakeFieldProcess.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchUtilities.h"

#include "AcceleratorModel/StdComponent/SWRFStructure.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "AcceleratorModel/StdComponent/SectorBend.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/utils.h"
#include "IO/MerlinIO.h"
#include "TLAS/TLASimp.h"

#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>


namespace {

#define COUT(x) cout<<std::setw(12)<<scientific<<std::setprecision(4)<<(x);

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;
using namespace ParticleTracking;

// needed to resolve gcc 3.2 ambiguity problem
inline double pow(int x, int y) { return pow(double(x),double(y)); }

Point2D GetSliceCentroid(ParticleBunch::const_iterator first,
                         ParticleBunch::const_iterator last) {
    Point2D c(0,0);
    double n=0;
    while(first!=last) {
        c.x += first->x();
        c.y += first->y();
        first++;
        n++;
    }
    return n>1 ? c/n : c;
}
PSvector GetSliceCentroid6D(ParticleBunch::const_iterator first,
                            ParticleBunch::const_iterator last) {
    PSvector c(0);
    double n=0;
    while(first!=last) {
        c += *first;
        first++;
        n++;
    }
    if(n>1)
        c/=n;
    return c;
}

};

namespace ParticleTracking {

CouplerWakeFieldProcess::CouplerWakeFieldProcess (int prio, size_t nb, double ns)
: WakeFieldProcess(prio,nb,ns,"COUPLERWAKEFIELD") {};
        
void CouplerWakeFieldProcess::SetCurrentComponent (AcceleratorComponent& component) {

    CombinedWakeRF* wake = dynamic_cast<CombinedWakeRF*>(component.GetWakePotentials());
    if(currentBunch!=0 && wake!=0) {
        clen = component.GetLength();
        switch(imploc) {
        case atCentre:
            impulse_s = clen/2.0;
            break;
        case atExit:
            impulse_s = clen;
            break;
        }
        current_s = 0;
        active = true;
        if(recalc || wake!=currentWake) {
            currentWake = wake;
            WakeFieldProcess::currentWake = currentWake;
            Init();
        }
    	const TWRFStructure* lcav = dynamic_cast< TWRFStructure* >(&component);
	if(lcav) {
	        phi=lcav->GetPhase();
        	V=lcav->GetVoltage();
        	//V=lcav->GetAmplitude();//Voltage/length
        	k=lcav->GetK();
	} else 	phi=V=k=0;
    }
    else {
        active = false;
        // check if we have a sector bend. If so, then
        // the bunch length will change and we need to rebin
        if(dynamic_cast<SectorBend*>(&component))
            recalc = true;
    }
}

void CouplerWakeFieldProcess::CalculateWakeT() {
// This routine is rather cpu intensive.
// First, calculate the transverse centroid of
// each bunch slice by taking the mean of the
// particle positions

    vector<Point2D> xyc;
    xyc.reserve(nbins+1);
    size_t i;
    for(i=0; i<nbins; i++)
        xyc.push_back(GetSliceCentroid(bunchSlices[i],bunchSlices[i+1]));
    xyc.push_back(xyc.back());
    // Now estimate the transverse bunch wake at the slice
    // boundaries in the same way we did for the longitudinal wake.

    double a0 = dz*(fabs(currentBunch->GetTotalCharge()))*ElectronCharge*Volt;
    wake_x = vector<double>(bunchSlices.size(),0.0);
    wake_y = vector<double>(bunchSlices.size(),0.0);

    //   a0=dz*Qtot*ElectronCharge*Volt
    //   volt=1.0e-9 (std unit is GeV) ElectronCharge=1.60219e-19C
    //   wake field unit V/C
    //   zmin = -nsig*sigz+z0;zmax =  nsig*sigz+z0;dz = (zmax-zmin)/nbins;

   for(i=0; i<bunchSlices.size(); i++) {
	
	   Vector2D cxy(0,0);
       if(i<bunchSlices.size()-1){
          // coupler wake kick (const, independent of z) at x[m],y[m]
          // cxy [V/m]
          cxy  = currentWake->Wxy(xyc[i].x,xyc[i].y);

	      // RF kick independent of bunch charge
	      // x[m],y[m],V[m]=Vacc
	      // phase[rad]=phi+k*dz,  dz[m]
	      // rfxy [1/cavity], V [GeV]
	      Vector2D rfxy 
		  = currentWake->CouplerRFKick(xyc[i].x,xyc[i].y,fabs(phi)+k*(zmin+(i+0.5)*dz));
//		  = currentWake->CouplerRFKick(xyc[i].x,xyc[i].y,fabs(phi));
	      wake_x[i] += rfxy.x*V/clen/a0;  // V[GeV] -> V[V]
	      wake_y[i] += rfxy.y*V/clen/a0;
       }
       for(int j=i; j<bunchSlices.size()-1; j++) {

	      // cavity transverse wake
	      // kick goes into the same direction as offset
	      // dx>0 => dx'>0
          double wxy = Qd[j]*(currentWake->Wtrans((j-i+0.5)*dz));
	      wake_x[i] += wxy*xyc[j].x;
	      wake_y[i] += wxy*xyc[j].y;

          // coupler kick assumed to be constant (short bunch) 
	      wake_x[i] += Qd[j]*cxy.x/clen;
	      wake_y[i] += Qd[j]*cxy.y/clen;  
	    
        }
        wake_x[i]*=a0;
        wake_y[i]*=a0;
    }

}


}; // end namespace ParticleTracking

