#include "BeamDynamics/CommonUtilities/BunchConverter.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamDynamics/SMPTracking/SliceMacroParticle.h"
#include "Random/RandomNG.h"
#include "NumericalUtils/utils.h"
#include <vector>

using namespace SMPTracking;
using namespace ParticleTracking;

// given an ParticleBunch we construct a SMPBunch with n_ct x n_dp particles
// adjust==true: force means(ParticleBunch)=means(SMPBunch) 
// nSigZ,nSigDP number of sigmas to scann
// nSigZ==0;nSigDP==0 scan from smallest to highest value in PB
// in both cases Qtot(PB)=Qtot(SMP)
// DK 1.4.2006
//
SMPTracking::SMPBunch* ParticleBunchConverter(ParticleTracking::ParticleBunch*  PB, int n_ct, int n_dp,
                                              double nSigZ, double nSigDP, bool adjust){

    double beamenergy = PB->GetReferenceMomentum();
    double Qtot       = PB->GetTotalCharge();
    PSmoments PSM;
    PB->GetMoments(PSM);// x xp y yp ct dp/p0
    double sig_dp = PSM.std(5);
    double sig_z  = PSM.std(4);

    // min,max is either min,max(PB) or n*sigmma
    double min_ct(0), max_ct(0);
    double min_dp(0), max_dp(0);
    if(nSigZ==0||nSigDP==0)
       for(ParticleBunch::const_iterator ipb = PB->begin(); ipb!= PB->end(); ipb++){
          double ct=ipb->ct();
          double dp=ipb->dp();
          if(ct<min_ct) min_ct = ct;
          if(ct>max_ct) max_ct = ct;
          if(dp<min_dp) min_dp = dp;
          if(dp>max_dp) max_dp = dp;
       }
    if(nSigZ!=0)  {max_ct = nSigZ*sig_z;    min_ct=-max_ct;}
    if(nSigDP!=0) {max_dp = nSigDP*sig_dp;  min_dp=-max_dp;}
 
    double dct=(max_ct-min_ct)/n_ct;
    double ddp=(max_dp-min_dp)/n_dp;
 
    // a histogram to integrate
    vector< vector<int> > h(n_ct, vector<int>(n_dp,0));
    double scale_ct=(n_ct-1)/(max_ct-min_ct);
    double scale_dp=(n_dp-1)/(max_dp-min_dp);
    for(ParticleBunch::const_iterator ipb = PB->begin(); ipb!= PB->end(); ipb++){
          int bin_ct=(ipb->ct()-min_ct)*scale_ct;
          int bin_dp=(ipb->dp()-min_dp)*scale_dp;
          // we assume a constant charge for all particles in a ParticleBunch !!
          if(bin_ct>=0&&bin_dp>=0&&bin_ct<n_ct&&bin_dp<n_dp) h[bin_ct][bin_dp]++;
    }
       
    SMPTracking::SMPBunch* bunch = new SMPTracking::SMPBunch(beamenergy,Qtot);

    double q0 = Qtot/PB->size();

    double ct  = min_ct+dct/2;
    for(int i_ct = 0; i_ct<n_ct; i_ct++, ct+=dct) {
       double dp  = min_dp+ddp/2;
       for(int i_dp = 0; i_dp<n_dp; i_dp++, dp+=ddp) {
           double q = q0*h[i_ct][i_dp];
           SliceMacroParticle m(PSM,ct,dp,q);
           bunch->AddParticle(m);
       }
    }
    
    return bunch;

};

