// Modified by D.Kruecker 18.2.2008
// to be used as base class for other wakefield types (spoiler,coupler,...)
//
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
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
                         ParticleBunch::const_iterator last)
{
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
                            ParticleBunch::const_iterator last)
{
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

WakeFieldProcess::WakeFieldProcess (int prio, size_t nb, double ns, string aID)
        : ParticleBunchProcess(aID,prio),imploc(atExit),nbins(nb),nsig(ns),currentWake(0),
        wake_x(0),wake_y(0),wake_z(0),Qd(),Qdp(),filter(0),recalc(true),inc_tw(true),
		oldBunchLen(0)
{
    SetFilter(14,2,1);
}

WakeFieldProcess::~WakeFieldProcess()
{
    if(filter)
        delete filter;
}

size_t WakeFieldProcess::CalculateQdist()
{
    pair<double,double> v = currentBunch->GetMoments(ps_CT);
    double z0 = v.first;
    double sigz = v.second;

    // calculate binning ranges
    zmin = -nsig*sigz+z0;
    zmax =  nsig*sigz+z0;
    dz = (zmax-zmin)/nbins;

    bunchSlices.clear();
    Qd.clear();
    Qdp.clear();

    // Qdp contains the slope of the charge distribution, smoothed using a filter
    size_t lost = ParticleBinList(*currentBunch,zmin,zmax,nbins,bunchSlices,Qd,Qdp,filter);

#ifndef NDEBUG
    ofstream os("qdist.dat");
    os<<zmin<<' '<<zmax<<' '<<dz<<endl;
    copy(Qd.begin(),Qd.end(),ostream_iterator<double>(os,"\n"));
#endif

    return lost;
}

// Smoothing filter takes the form of a set of coefficients
// calculated using the Savitzky-Golay technique
// n gives the width of the window on either side of the reference point
// m gives the order of the polynomial fitted to the points within the window
// d gives the order of the derivative required
// For CSR wake we need the first derivative
void WakeFieldProcess::SetFilter(int n, int m, int d)
{
    if(filter)
        delete filter;

    filter = new vector<double>;
    savgol(*filter, n, n, d, m);

#ifndef NDEBUG
    ofstream os("filter.dat");
    copy(filter->begin(),filter->end(),ostream_iterator<double>(os,"\n"));
#endif
}

void WakeFieldProcess::SetCurrentComponent (AcceleratorComponent& component)
{
    //	TWRFStructure* cavity = dynamic_cast<TWRFStructure*>(&component);
    //	WakePotentials* wake = cavity!=0 ? cavity->GetWakePotentials() : 0;


    WakePotentials* wake = component.GetWakePotentials();
    // if not initialize(=0) we assume that
    // WakeFieldProcess is responsible - for backward compatibility
    // in general expected process must be equal to this process
    if( wake &&
        wake->GetExpectedProcess()!=0 &&
        typeid(*(wake->GetExpectedProcess()))!=typeid(*this))  
        wake=0;

    //if(wake!=0) cout<<GetID()<<endl;

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
    // check if we are responsible for the current wakefield
    // a class derived from this class must 
    // include the appropriate check for its own wake potential type!
    if(typeid(WakePotentials*)!=typeid(currentWake)) return;
    
    // here we apply the wake field for
    // the step ds
    size_t n=0;
    double p0 = currentBunch->GetReferenceMomentum();

    // If the bunch length or binning has been changed,
    // we must recalculate the wakes
	// dk explicit check on bunch length
	if(recalc||oldBunchLen!=currentBunch->size())
        Init();

    // We always need to recalculate the transverse wake
    if(inc_tw)
        CalculateWakeT();

    // Now iterate over the bunch slices,
    // calculating the wakefield kicks using
    // linear interpolation between the values
    // at the slice boundaries

    double bload=0;

#define WAKE_GRADIENT(wake) (wake).empty() ? 0 : ((wake[nslice+1]-wake[nslice])/dz);

    double z=zmin;

    for(size_t nslice = 0; nslice<nbins; nslice++) {

        double gz = WAKE_GRADIENT(wake_z);
        double gx = WAKE_GRADIENT(wake_x);
        double gy = WAKE_GRADIENT(wake_y);


 

       for(ParticleBunch::iterator p=bunchSlices[nslice]; p!=bunchSlices[nslice+1]; p++) {
            double zz = p->ct()-z;
            double ddp = -ds*(wake_z[nslice]+gz*zz)/p0;
            p->dp() += ddp;
            bload += ddp;

            double dxp =  inc_tw? ds*(wake_x[nslice]+gx*zz)/p0 : 0;
            double dyp =  inc_tw? ds*(wake_y[nslice]+gy*zz)/p0 : 0;
           

            p->xp() = (p->xp()+dxp)/(1+ddp);
            p->yp() = (p->yp()+dyp)/(1+ddp);
        }
        z+=dz;
    }
    if(!currentWake->Is_CSR())
        currentBunch->AdjustRefMomentum(bload/currentBunch->size());
    
}

double WakeFieldProcess::GetMaxAllowedStepSize () const
{
    return impulse_s-current_s;
}

void WakeFieldProcess::Init()
{
    double Qt  = currentBunch->GetTotalCharge();
	
	//keep track of bunch length to be aware of modifications
	oldBunchLen=currentBunch->size();

    size_t nloss = CalculateQdist();
    if(nloss!=0) {
        // Even though we have truncated particles, we still keep the
        // the bunch charge constant
        currentBunch->SetMacroParticleCharge(Qt/(currentBunch->size()));
        MerlinIO::warning()<<GetID()<<" (WakefieldProcess): "<<nloss<<" particles truncated"<<endl;
    }

    // Calculate the long. bunch wake.
    CalculateWakeL();
    recalc=false;
}



void WakeFieldProcess::CalculateWakeL()
{
    wake_z = vector<double>(bunchSlices.size(),0.0);
    double a0 = dz*fabs(currentBunch->GetTotalCharge())*ElectronCharge*Volt;

    // Estimate the bunch wake at the slice boundaries by
    // convolving the point-like wake over the current bunch
    // charge distribution Qd.
    //
    // Note that the distribution Qd is estimated at the
    // centre of each slice, not the slice boundary.
    //
    // CSR wake differs from classical wakes in two respects:
    // 1) wake from the tail of the bunch affects the head
    // 2) amplitude of wake depends on slope of charge distribution
    //    rather than directly on the distribution
    // Code to handle CSR wake added by A.Wolski 12/2/2003

    if(currentWake->Is_CSR())
    {
        for(int i=0; i<bunchSlices.size(); i++) {
            for(int j=1; j<i; j++) {
                wake_z[i] += Qdp[j]*(currentWake->Wlong((j-i+0.5)*dz))/dz;
            }
            wake_z[i]*=a0;
        }
    }
    else
    {
        for(int i=0; i<bunchSlices.size(); i++) {
            for(int j=i; j<bunchSlices.size()-1; j++) {
                wake_z[i] += Qd[j]*(currentWake->Wlong((j-i+0.5)*dz));
            }
            wake_z[i]*=a0;
        }
    }

#ifndef NDEBUG
    ofstream os("bunchWake.dat");
    os<<zmin<<' '<<zmax<<' '<<dz<<endl;
    copy(wake_z.begin(),wake_z.end(),ostream_iterator<double>(os,"\n"));
#endif

}

void WakeFieldProcess::CalculateWakeT()
{
    // This routine is rather cpu intensive.
    // First, calculate the transverse centroid of
    // each bunch slice by taking the mean of the
    // particle positions

    vector<Point2D> xyc;
    xyc.reserve(nbins);
    size_t i;
    for(i=0; i<nbins; i++)
        xyc.push_back(GetSliceCentroid(bunchSlices[i],bunchSlices[i+1]));

    // Now estimate the transverse bunch wake at the slice
    // boundaries in the same way we did for the longitudinal wake.

    double a0 = dz*(fabs(currentBunch->GetTotalCharge()))*ElectronCharge*Volt;
    wake_x = vector<double>(bunchSlices.size(),0.0);
    wake_y = vector<double>(bunchSlices.size(),0.0);
    for(i=0; i<bunchSlices.size(); i++) {
        for(int j=i; j<bunchSlices.size()-1; j++) {
            double wxy = Qd[j]*(currentWake->Wtrans((j-i+0.5)*dz));
            wake_x[i] += wxy*xyc[j].x;
            wake_y[i] += wxy*xyc[j].y;
        }
        wake_x[i]*=a0;
        wake_y[i]*=a0;
    }
}

void WakeFieldProcess::DumpSliceCentroids(ostream& os) const
{
    for(size_t i=0; i<nbins; i++) {
        os<<std::setw(4)<<i;
        os<<GetSliceCentroid6D(bunchSlices[i],bunchSlices[i+1]);
    }
}

void WakeFieldProcess::InitialiseProcess (Bunch& bunch)
{
    ParticleBunchProcess::InitialiseProcess(bunch);
    currentWake = 0;
    recalc = true;
}

// Calculate the Savitsky-Golay smoothing filter
// Adapted from Numerical Recipes in C
// This routine need only be executed once,
// when the WakeFieldProcess is initialized.
int powi(int i, int j)
{
	int p = 1;
	for(int m=0; m<j; m++)
		p *= i;
	return p;
};

void savgol(vector<double>& c, int nl, int nr, int ld, int m)
{
    Matrix<double> a(m+1,m+1);

    for(int i=0; i<=m; i++)
        for(int j=0; j<=i; j++)
        {
            double sum = 0.0;

            for(int k=-nl; k<=nr; k++)
                sum += powi(k, i) * powi(k, j);

            a(i,j) = sum;
            a(j,i) = sum;
        }

    vector<int> indx;
    double d;
    ludcmp(a,indx,d);

    Vector<double> b(m+1);
    for(int j=0; j<=m; j++)
        b(j) = 0.0;
    b(ld) = 1.0;

    lubksb(a,indx,b);

    c.clear();

    for(int k=-nl; k<=nr; k++)
    {
        double sum = b(0);
        double fac = 1.0;

        for(int mm=1; mm<=m; mm++)
            sum += b(mm) * (fac *= k);

        c.push_back(sum);
    }
}

}; // end namespace ParticleTracking

