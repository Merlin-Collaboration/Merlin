/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/01 12:36:35 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#include "BeamDynamics/SMPTracking/SMPBunch.h"
#include "BeamDynamics/SMPTracking/SMPTransform3D.h"
#include <fstream>
#include <iterator>

namespace SMPTracking {

using namespace std;

SMPBunch::SMPBunch (double p, double q)
        : Bunch(p,q),Qt(0)
{}

SMPBunch::SMPBunch (const std::string& fname)
        : Bunch(1,1),Qt(0)
{
    ifstream ifs(fname.c_str());
    if(!ifs) {
        cerr<<"SMPBunch file not found: "<<fname<<endl;
        abort();
    }

    double p0,ct0;
	size_t nmp;
    if(ifs>>p0>>ct0>>nmp) {
        SetReferenceMomentum(p0);
        SetReferenceTime(ct0);
        slices.reserve(nmp);
    }
    else {
        cerr<<"Bad SMPBunch file format: "<<fname<<endl;
        abort();
    }
    while(ifs) {
        SliceMacroParticle p;
        p.Read(ifs);
        if(ifs) {
            slices.push_back(p);
            Qt+=p.Q();
            nmp--;
        }
    }

    if(nmp!=0) {
        cerr<<"Bad SMPBunch file format: "<<fname<<endl;
        abort();
    }

    sort(begin(),end());
    SetChargeSign(Qt);
};

SMPBunch::~SMPBunch ()
{
    // nothing to do
}

double SMPBunch::GetTotalCharge () const
{
    return Qt;
}

PSmoments& SMPBunch::GetMoments (PSmoments& sigma) const
{
    sigma.zero();
    int i,j;
    for(const_iterator p=begin(); p!=end(); p++) {
        const SliceMacroParticle& x = (*p);
        double w = x.Q()/Qt;
        for(i=0; i<6; i++) {
            sigma[i]+=w*x[i];  // centroid
            for(j=0; j<=i; j++) //2nd-order moments
                sigma(i,j)+= w*(x[i]*x[j]+x(i,j));
        }
    }
    for(i=0; i<6; i++)
        for(j=0; j<=i; j++)
            sigma(i,j)-=sigma[i]*sigma[j];

    return sigma;
}

PSmoments2D& SMPBunch::GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const
{
    for(const_iterator p=begin(); p!=end(); p++) {
        const SliceMacroParticle& x = (*p);
        double w = x.Q()/Qt;

        sigma[0]+=w*x[u];
        sigma[1]+=w*x[v];

        sigma(0,0)+=w*(x(u,u)+x[u]*x[u]);
        sigma(0,1)+=w*(x(u,v)+x[u]*x[v]);
        sigma(1,1)+=w*(x(v,v)+x[v]*x[v]);
    }
    sigma(0,0)-=sigma[0]*sigma[0];
    sigma(0,1)-=sigma[0]*sigma[1];
    sigma(1,1)-=sigma[1]*sigma[1];

    return sigma;

}

PSvector& SMPBunch::GetCentroid (PSvector& x) const
{
    x.zero();
    for(const_iterator p=begin(); p!=end(); p++)
        x+=p->GetChargeWeightedCentroid();
    x/=Qt;
    return x;
}

Point2D SMPBunch::GetProjectedCentroid (PScoord u, PScoord v) const
{
    Point2D x0(0,0);
    for(const_iterator p=begin(); p!=end(); p++)
        x0+=p->GetChargeWeightedCentroid(u,v);
    x0/=Qt;
    return x0;
}

double SMPBunch::AdjustRefMomentumToMean ()
{
    double dp0=0;
    for(const_iterator cp=begin(); cp!=end(); cp++)
        dp0+=(cp->dp())*(cp->Q());
    double ddp=dp0/Qt+1.0;
    for(iterator p=begin(); p!=end(); p++)
        p->dp()/=ddp;
    SetReferenceMomentum(GetReferenceMomentum()*ddp);
    return GetReferenceMomentum();
}

double SMPBunch::AdjustRefTimeToMean ()
{
    double ct=0;
    for(const_iterator p=begin(); p!=end(); p++)
        ct+=(p->ct())*(p->Q());
    IncrReferenceTime(ct/Qt);
    return GetReferenceTime();
}

void SMPBunch::Output (std::ostream& os) const
{
    static char delim[]=" ";

    os<<GetReferenceMomentum()<<delim;
    os<<GetReferenceTime()<<delim;
    os<<slices.size()<<endl;
    copy(begin(),end(),ostream_iterator<SliceMacroParticle>(os,""));
}

Histogram& SMPBunch::ProjectDistribution (PScoord axis, Histogram& hist) const
{
    // TO DO
    return hist;
}

bool SMPBunch::ApplyTransformation (const Transform3D& t)
{
    if(t.isIdentity())
        return true;

    SMPTransform3D xf(t);
    for(iterator p=begin(); p!=end(); p++)
        xf.Apply(*p);

    return false;
}

void SMPBunch::AddParticle(const SliceMacroParticle& p, bool do_sort)
{
    slices.push_back(p);
    if(do_sort)
        sort(begin(),end());
    Qt+=p.Q();
}

void SMPBunch::SortByCT()
{
    sort(begin(),end());
}

void SMPBunch::AdjustCentroid(const PSvector & s)
{
    for(iterator p=begin(); p!=end(); p++) {
        (*p).x() += s.x();
        (*p).xp() += s.xp();
        (*p).y() += s.y();
        (*p).yp() += s.yp();
        (*p).ct() += s.ct();
        (*p).dp() += s.dp();
    }
}

}; // end namespace SMPTracking
