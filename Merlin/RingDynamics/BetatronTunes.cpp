/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 23:03:54 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/CollimateParticleProcess.h"
#include "NumericalUtils/NumericalConstants.h"

#include "BetatronTunes.h"

using namespace std;
using namespace ParticleTracking;

BetatronTunes::BetatronTunes(AcceleratorModel* aModel, double refMomentum)
        : theModel(aModel), p0(refMomentum) {}

void BetatronTunes::FindTunes(PSvector& particle, int ntrack, bool diffusion)
{
    vector<double> xData1;
    xData1.reserve(ntrack * 2);

    vector<double> yData1;
    yData1.reserve(ntrack * 2);

    vector<double> xData2;
    xData2.reserve(ntrack * 2);

    vector<double> yData2;
    yData2.reserve(ntrack * 2);

    ParticleBunch bunch(p0,1.0);
    bunch.push_back(particle);

    ParticleTracker tracker(theModel->GetBeamline(), &bunch, false);

	CollimateParticleProcess collimate(1,COLL_AT_EXIT);
	collimate.SetLossThreshold(200);
	tracker.AddProcess(&collimate);

    vector<double>* xData = &xData1;
    vector<double>* yData = &yData1;

    bool stable = true;
    int trackn = diffusion ? ntrack*2 : ntrack;
    for(int nturn=0; nturn<trackn; nturn++)
    {
        if(nturn==0)
            tracker.Run();
        else
            tracker.Continue();
        ParticleBunch& tracked_bunch = tracker.GetTrackedBunch();
		if(tracked_bunch.size()==0){
			stable = false;
			break;
		}
        ParticleBunch::iterator pp = tracked_bunch.begin();
        Particle& p1 = *pp;
#ifdef _MSV_VER
        if( _isnan(p1.x()) || _isnan(p1.y()) )
#else
        if( fabs(p1.x())>1.0e+10 || fabs(p1.y())>1.0e+10 )
#endif
        {
            stable = false;
            break;
        }

        double hanningFilter = sin( pi*double(nturn+1)/double(ntrack) );
        hanningFilter *= hanningFilter;

        if(nturn==ntrack)
        {
            xData = &xData2;
            yData = &yData2;
        }

        xData->push_back( p1.x()  * hanningFilter );
        xData->push_back( p1.xp() * hanningFilter );
        yData->push_back( p1.y()  * hanningFilter );
        yData->push_back( p1.yp() * hanningFilter );
    }

    Qx = stable ? FindTune(xData1) : 0;
    Qy = stable ? FindTune(yData1) : 0;

    double Qx2 = (diffusion & stable) ? FindTune(xData2) : 0;
    double Qy2 = (diffusion & stable) ? FindTune(yData2) : 0;

    dQx = Qx2 - Qx;
    dQy = Qy2 - Qy;
}

double BetatronTunes::FindTune(vector<double>& data)
{
    vector<double> spectrum;
    int nvals = data.size()/2;
    spectrum.reserve(nvals);

    FFT(data);
    double peak_height = 0.0;
    int peak_position = 0;

    for(int i=0; i<2*nvals; i+=2)
    {
        double vamp = sqrt(data[i]*data[i] + data[i+1]*data[i+1]);
        spectrum.push_back(vamp);
        if(vamp > peak_height)
        {
            peak_height = vamp;
            peak_position = i/2;
        }
    }

    double norm = sqrt(double(nvals));
    double centr_peak = spectrum[peak_position]/norm;
    double right_peak = spectrum[peak_position+1]/norm;

    double peak_correction =
        nvals
        * asin( amp( centr_peak, right_peak, cos(twoPi/nvals) ) * sin(twoPi/double(nvals)) )
        / twoPi;

    return (peak_position + peak_correction)/nvals;
}

void BetatronTunes::FFT(vector<double>& data)
{
    int nn=data.size()/2;
    int n = nn << 1;
    int j = 1;
    for(int i=1; i<n; i+=2)
    {
        if(j>i)
        {
            swap(data[j-1],data[i-1]);
            swap(data[j],data[i]);
        }

        int m = nn;
        while(m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    int mmax = 2;
    while(n > mmax)
    {
        int istep = mmax << 1;
        double theta = twoPi / mmax;
        double wtemp = sin(theta/2.);
        double wpr = -2.0 * wtemp * wtemp;
        double wpi = sin(theta);
        double wr = 1.;
        double wi = 0.;
        for(int m=1; m<mmax; m+=2)
        {
            for(int i=m; i<=n; i+=istep)
            {
                j = i + mmax;
                double tempr = wr * data[j-1] - wi * data[j];
                double tempi = wr * data[j] + wi * data[j-1];
                data[j-1] = data[i-1] - tempr;
                data[j] = data[i] - tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wr = (wtemp=wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
    }

}

double BetatronTunes::amp(double a, double b, double c)
{
    double a1 = -(a+b*c)*(a-b);
    double b1 = b * sqrt( c*c*(a+b)*(a+b) - 2.*a*b*(2*c*c-c-1.) );
    double c1 = a*a + b*b + 2.*a*b*c;
    return (a1+b1)/c1;
}
