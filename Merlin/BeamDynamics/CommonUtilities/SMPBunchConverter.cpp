#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamModel/PSvector.h"
#include "BeamDynamics/CommonUtilities/BunchConverter.h"
#include <iostream>
#include "NumericalUtils/TCovMtrx.h"
#include "Random/MultiNormal.h"
#include <set>

using namespace ParticleTracking;
using namespace SMPTracking;

// Given an SPMBunch we construct a ParticleBunch with about N particles.  
// In general the final number of particles is !=N (equal N only on average).
// For adjust=true (default) the centroid (6d) of the particle bunch is adjusted to the 
// SMP bunch centroid.
// The discrete ct, dp points in the original SMPBunch are smeared out in the ParticleBunch 
// with gauss(delta/2) where delta is given by the distance in either ct or dp.
// DK 1.4.2006
//
ParticleTracking::ParticleBunch* SMPBunchConverter(SMPTracking::SMPBunch*  SB, int N,bool adjust){

    if(N<SB->Size()) {
       cout<<"SMPBunchConverter: Warning, you asked for a low number of particles!"<<endl; 
       //return new ParticleBunch(0);
    }
    double beamenergy=SB->GetReferenceMomentum();
    double reftime=SB->GetReferenceTime();
    double Qtot = SB->GetTotalCharge();

    // to calculate distances in ct and dp in the main loop
    set<double> set_ct;
    set<double> set_dp;
    for(SMPBunch::const_iterator sp=SB->begin(); sp!=SB->end(); sp++) {
           set_ct.insert(sp->ct());
           set_dp.insert(sp->dp());
    }
    int N_ct=set_ct.size();
    int N_dp=set_dp.size();

    PSvector mean(0);
    PSvectorArray particles;
    PSvector p;
    for(SMPBunch::const_iterator sp=SB->begin(); sp!=SB->end(); sp++) {
        const SliceMacroParticle& x = (*sp);
        //x.Write(cout);

        // the 4 dimensinal normal random generator
        MultiNormal<4> MN(x);

        // weight for this slice
        double w = x.Q()/Qtot;

        // # of particles to represent this slice
        // sum_i(int(N*w_i))<=N, <sum_i(int(N*w_i))>=N-k/2 , sum_i(w_i)=1; i=1..k
        int    n = N*w+0.5;    //  <sum_i(n)> = N

        double sct,sdp; // sigmas for smearing - sig/2 for a smooth distrib
        // 1/2(x_n-x_(n+1) or x_n-x_(n-1))
        set<double>::iterator itPlus;
        if(N_ct>1){ 
            itPlus = set_ct.upper_bound(x.ct());
            if(itPlus==set_ct.end()) sct = fabs(*--set_ct.rbegin()-x.ct())/2;
            else                     sct = fabs(*itPlus-x.ct())/2;
        } else sct=0;
        sct*=sct;
        if(N_dp>1){ 
            itPlus = set_dp.upper_bound(x.dp());
            if(itPlus==set_dp.end()) sdp = fabs(*--set_dp.rbegin()-x.dp())/2;
            else                     sdp = fabs(*itPlus-x.dp())/2;
        } else sdp=0;
        sdp*=sdp;

        for(int i=0;i<n;i++){
            RealVector v = MN.GetRandVec();
            p.x()  = v(0);
            p.xp() = v(1);
            p.y()  = v(2);
            p.yp() = v(3);
            p.ct() = sct>0?RandomNG::normal(x.ct(),sct):x.ct();
            p.dp() = sdp>0?RandomNG::normal(x.dp(),sdp):x.dp();
            mean+=p;
            particles.push_back(p);
        }
    }

    if(adjust){
       mean/=particles.size();
       PSvector delta = SB->GetCentroid(delta)-mean;
       for(PSvectorArray::iterator it=particles.begin();it!=particles.end();it++) *it+=delta;
    }

    ParticleBunch* np = new ParticleBunch (beamenergy, Qtot, particles);
    np->SetReferenceTime(reftime);

    return np;

};
