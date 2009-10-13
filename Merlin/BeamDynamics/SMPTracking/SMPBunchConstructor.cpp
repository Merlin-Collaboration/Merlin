/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/08/02 15:16:18 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#include "BeamDynamics/SMPTracking/SMPBunchConstructor.h"
#include "BasicTransport/NormalTransform.h"
#include "NumericalUtils/utils.h"

namespace {

void MakeDistribution(double ns, vector<double>& d)
{
    double dx = 2*ns/d.size();
    double x = -ns;
    for(size_t n=0; n<d.size(); n++, x+=dx)
        d[n] = NormalBin(x,x+dx);
}

}; // end anonymous namespace


namespace SMPTracking {


SMPBunchConstructor::SMPBunchConstructor (const BeamData& beam, size_t ns1, size_t nsm)
        : ns(ns1),np(nsm),beamdat(beam),nSigZ(3.0),nSigDP(3.0)
{}

SMPBunchConstructor::~SMPBunchConstructor ()
{}

Bunch* SMPBunchConstructor::ConstructBunch (int) const
{
    SMPBunch* bunch = new SMPBunch(beamdat.p0,beamdat.charge);

    // First, store the dispersion correlations and set them to
    // zero.
    BeamData bd = beamdat;
    bd.Dx=bd.Dxp=bd.Dy=bd.Dyp=0.0;
    double dx = beamdat.Dx;
    double dxp = beamdat.Dxp;
    double dy = beamdat.Dy;
    double dyp = beamdat.Dyp;
    PSmoments S = BeamDataToSigmaMtrx(bd);

    // the macroparticles are constructed to represent the
    // total charge in a bin of width 6sigma/ns

    double q0 = beamdat.charge;
    double qt = 0;

    vector<double> zdist(ns);
    vector<double> dpdist(np);

    MakeDistribution(nSigZ,zdist);
    MakeDistribution(nSigDP,dpdist);

    double dz  = 2*nSigZ*beamdat.sig_z/ns;
    double ddp = 2*nSigDP*beamdat.sig_dp/np;

    double z   =-nSigZ*beamdat.sig_z+dz/2;
    for(size_t nSlice = 0; nSlice<ns; nSlice++, z+=dz) {
        double dp = -nSigDP*beamdat.sig_dp+ddp/2.0;
        for(size_t nPart = 0; nPart<np; nPart++,dp+=ddp) {
            double q = q0*zdist[nSlice]*dpdist[nPart];
            SliceMacroParticle m(S,z,dp,q);
            // Add dispersion correlations to centroid
            m.x() += dx*dp;
            m.xp()+= dxp*dp;
            m.y() += dy*dp;
            m.yp()+= dyp*dp;
            bunch->AddParticle(m);
            qt+=q;
        }
    }
    //cout<<"total generated charge = "<<qt<<endl;
    PSvector x0;
    bunch->GetCentroid(x0);
    return bunch;
}

SMPBunch* SMPBunchConstructor::ConstructSMPBunch () const
{
    return static_cast<SMPBunch*>(ConstructBunch());
}

}; // end namespace SMPTracking

