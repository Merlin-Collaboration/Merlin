// WakeFieldProcess implementation
//
#include "BeamDynamics/SMPTracking/SMPWakeFieldProcess.h"
#include "AcceleratorModel/WakePotentials.h"
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

// unit conversion factor for wake potential
using namespace PhysicalConstants;
using namespace PhysicalUnits;
using namespace SMPTracking;
const double wconv = ElectronCharge*Volt;


Point2D GetSliceCentroid(SMPBunch::const_iterator first, SMPBunch::const_iterator last)
{
    double qt=0;
    Point2D X(0,0);
    while(first!=last) {
        double q = first->Q();
        X.x += (first->x())*q;
        X.y += (first->y())*q;
        qt+=q;
        first++;
    }
    X.x/=qt;
    X.y/=qt;
    return X;
}

// Apply Wake Kick (AWK). Used to apply the same kick to a range of macro-particles
struct AWK {

    double Eref_0;
    double Eref_1;
    double dPx;
    double dPy;
    double dPz;

    AWK(double E0, double dE, double px, double py, double pz)
            : Eref_0(E0),Eref_1(E0+dE),dPx(px),dPy(py),dPz(pz) {}

    void operator()(SliceMacroParticle& p) const {
        double Ei=Eref_0*(1+p.dp());
        double Ef=Ei+dPz;
        double Er=Ei/Ef;
        double Er2=Er*Er;

        p.xp()=(Ei*p.xp()+dPx)/Ef;
        p.yp()=(Ei*p.yp()+dPy)/Ef;
        p.dp()=Ef/Eref_1-1;

        // Scale second-order moments
        p(0,1)*=Er;
        p(0,3)*=Er;
        p(1,1)*=Er2;
        p(1,3)*=Er2;
        p(2,3)*=Er;
        p(3,3)*=Er2;
    }
};
};

namespace SMPTracking {


WakeFieldProcess::WakeFieldProcess (int prio, double sw, string aID)
        : SMPBunchProcess(aID,prio),imploc(atExit),currentWake(0),recalc(true),inc_tw(true),dz(sw)
{}

WakeFieldProcess::~WakeFieldProcess()
{}


void WakeFieldProcess::SetCurrentComponent (AcceleratorComponent& component)
{
    WakePotentials* wake = component.GetWakePotentials();
    // if not initialize(=0) we assume that
    // WakeFieldProcess is responsible - for backward compatibility
    // in general expected process must be equal to this process
    if( wake &&
        wake->GetExpectedProcess()!=0 &&
        typeid(*(wake->GetExpectedProcess()))!=typeid(*this))
        wake=0;


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
            Init();
        }
    }
    else {
        active = false;
        // check if we have a sector bend. If so, then
        // the bunch length will change and we need to rebin
        if(dynamic_cast<SectorBend*>(&component))
            recalc = true;
    }
}

void WakeFieldProcess::DoProcess(double ds)
{
    current_s+=ds;
    if(fequal(current_s,impulse_s)) {
        ApplyWakefield(clen);
        active = false;
    }
}

void WakeFieldProcess::ApplyWakefield(double ds)
{
    // here we apply the wake field for
    // the step ds
    size_t n=0;
    double p0 = currentBunch->GetReferenceMomentum();

    // If the bunch length or binning has been changed,
    // we must recalculate the wakes
    if(recalc)
        Init();

    size_t np = slice_z.size();
    double dE= -bload*ds;

    for(size_t i=0; i<np; i++) {

        double kick_x = 0;
        double kick_y = 0;

        if(inc_tw) {
            for(size_t j=i+1; j<np; j++) {
                double w = (currentWake->Wtrans(slice_z[j]-slice_z[i]))*slice_q[j];
                Point2D X0=GetSliceCentroid(sliceBoundaries[j],sliceBoundaries[j+1]);
                kick_x += w*X0.x;
                kick_y += w*X0.y;
            }
        }

        double dPx =  ds*wconv*kick_x;
        double dPy =  ds*wconv*kick_y;
        double dPz = -ds*wake_z[i]; // beamloading

        // apply the kicks to the macro-particles in the i-th slice
        for_each(sliceBoundaries[i],sliceBoundaries[i+1],AWK(p0,dE,dPx,dPy,dPz));
    }
    currentBunch->SetReferenceMomentum(p0+dE);
}

double WakeFieldProcess::GetMaxAllowedStepSize () const
{
    return impulse_s-current_s;
}

void WakeFieldProcess::Init()
{
    PrepSlices();
    PrepLWake();
    recalc=false;
}

void WakeFieldProcess::PrepLWake()
{
    // Calculate and cache in wake_z the longitudinal bunch
    // wake per unit length.
    size_t np = slice_z.size();
    wake_z = vector<double>(np,0.0);
    double qt=0;
    bload=0;
    for(size_t i=0; i<np; i++) {
        for(size_t j=i; j<np; j++) {
            double dz = slice_z[j]-slice_z[i];
            // note factor 0.5 for FTBL
            wake_z[i]+=(i==j?0.5:1.0)*(currentWake->Wlong(dz))*slice_q[j]*wconv;
        }
        // beam loading
        bload+=slice_q[i]*wake_z[i];
        qt+=slice_q[i];
    }
    bload/=qt; // in GeV/m
    //cout<<"beam loading = "<<bload/keV<<" keV/m"<<endl;
}

void WakeFieldProcess::InitialiseProcess (Bunch& bunch)
{
    SMPBunchProcess::InitialiseProcess(bunch);
    currentWake = 0;
    recalc = true;
}

void WakeFieldProcess::PrepSlices()
{
    vector<SMPBunch::iterator> sb;
    vector<double> tz;
    vector<double> tq;
    sb.reserve(32);
    tz.reserve(32);
    tq.reserve(32);

    SMPBunch::iterator p = currentBunch->begin();
    sb.push_back(p);

    double z=p->ct()+dz;
    double zavg=(p->ct())*(p->Q());
    double qt=p->Q();

    for(;p!=currentBunch->end(); p++) {
        if(p->ct()>=z) {
            sb.push_back(p);
            tz.push_back(zavg/qt);
            tq.push_back(qt);
            z=p->ct()+dz;
            zavg=(p->ct())*(p->Q());
            qt=p->Q();
        }
        else {
            zavg+=(p->ct())*(p->Q());
            qt+=p->Q();
        }
    }

    sb.push_back(p); // end()
    tz.push_back(zavg/qt);
    tq.push_back(qt);

    sliceBoundaries.swap(sb);
    slice_z.swap(tz);
    slice_q.swap(tq);

    //		cout<<slice_z.size()<<" bins identified"<<endl;
}

}; // end namespace SMPTracking
