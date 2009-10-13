/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include <list>
#include <iterator>
// Transform3D
#include "EuclideanGeometry/Transform3D.h"
// PSvectorTransform3D
#include "BasicTransport/PSvectorTransform3D.h"
// ParticleBunch
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;

// needed by sort algorithm
inline bool operator<(const PSvector& p1, const PSvector& p2)
{
    return p1.ct()<p2.ct();
}

namespace {

using namespace ParticleTracking;

template<class V>
struct APSV1 {
    V& p;
    const double w;

    APSV1(V& p0, double n) : p(p0),w(1/n) {}
    void operator()(const PSvector& p1) {
        for(int i=0; i<6; i++) p[i]+=w*p1[i];
    }
};

struct APSV2 {
    const int u,v;
    Point2D X;
    const double w;

    APSV2(int i,int j,double n) : u(i),v(j),X(0,0),w(1/n){}
    void operator()(const PSvector& p) {
        X.x+=w*p[u];
        X.y+=w*p[v];
    }
};

struct PSVVAR1 {
    PSmoments& S;
    const double w;
    PSVVAR1(PSmoments& sig,double n):S(sig),w(1/n) {}
    void operator()(const PSvector& p) {
        for(int i=0; i<6; i++)
            for(int j=0; j<=i; j++)
                S(i,j)+=w*(p[i]-S[i])*(p[j]-S[j]);
    }
};

struct PSVVAR2 {
    PSmoments2D& S;
    const int u,v;
    const double w;

    PSVVAR2(int u1, int v1, double n, PSmoments2D& sig)
            :S(sig),u(u1),v(v1),w(1/n) {}

    void operator()(const PSvector& p) {
        double x=p[u]-S[0];
        double y=p[v]-S[1];
        S(0,0)+=w*x*x;
        S(0,1)+=w*x*y;
        S(1,1)+=w*y*y;
    }
};

template<class IT>
double Mean(IT F, IT L,int i)
{
    double s=0;
    double n=0;
    while(F!=L) {s+=(*F++)[i];n++;}
    return s/n;
}

template<class T>
inline void SortArray(std::vector<T>& array)
{
    sort(array.begin(),array.end());
}

template<class T>
inline void SortArray(std::list<T>& array)
{
    array.sort();
}

};

namespace ParticleTracking {

ParticleBunch::ParticleBunch (double P0, double Q, PSvectorArray& particles)
        : Bunch(P0,Q),qPerMP(Q/particles.size()),pArray()
{
    pArray.swap(particles);
}

ParticleBunch::ParticleBunch (double P0, double Q, std::istream& is)
        : Bunch(P0,Q)
{
    PSvector p;
    while(is>>p)
        push_back(p);

    qPerMP = Q/size();
}

ParticleBunch::ParticleBunch (double P0, double Qm)
        : Bunch(P0,Qm),qPerMP(Qm)
{}

double ParticleBunch::GetTotalCharge () const
{
    return qPerMP*size();
}

PSmoments& ParticleBunch::GetMoments (PSmoments& sigma) const
{
    sigma.zero();
    for_each(begin(),end(),APSV1<PSmoments>(sigma,size()));
    for_each(begin(),end(),PSVVAR1(sigma,size()));
    return sigma;
}

PSmoments2D& ParticleBunch::GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const
{
    Point2D X = for_each(begin(),end(),APSV2(u,v,size())).X;
    sigma[0]=X.x;
    sigma[1]=X.y;
    for_each(begin(),end(),PSVVAR2(u,v,size(),sigma));
    return sigma;
}

PSvector& ParticleBunch::GetCentroid (PSvector& p) const
{
    p.zero();
    for_each(begin(),end(),APSV1<PSvector>(p,size()));
    return p;
}

std::pair<double,double> ParticleBunch::GetMoments(PScoord i) const
{
    double u = Mean(begin(),end(),i);
    double v=0;
    for(const_iterator p = begin(); p!=end(); p++) {
        double x = (*p)[i]-u;
        v+=x*x;
    }

    return make_pair(u,sqrt(v/size()));
}

Point2D ParticleBunch::GetProjectedCentroid (PScoord u, PScoord v) const
{
    return for_each(begin(),end(),APSV2(u,v,size())).X;
}

double ParticleBunch::AdjustRefMomentumToMean ()
{
    return AdjustRefMomentum( Mean(begin(),end(),ps_DP) );
}

double ParticleBunch::AdjustRefMomentum (double dpp)
{
    double onePlusDpp = 1+dpp;

    for(iterator p=begin(); p!=end(); p++)
        p->dp() = (p->dp()-dpp)/onePlusDpp;

    double P0 = onePlusDpp*GetReferenceMomentum();
    SetReferenceMomentum(P0);
    return P0;
}

double ParticleBunch::AdjustRefTimeToMean ()
{
    double meanct = Mean(begin(),end(),ps_CT);
    for(iterator p=begin(); p!=end(); p++)
        (*p).ct()-=meanct;

    double CT = GetReferenceTime()-meanct;
    SetReferenceTime(CT);
    return CT;
}

Histogram& ParticleBunch::ProjectDistribution (PScoord axis, Histogram& hist) const
{
    // TODO:
    return hist;
}


bool ParticleBunch::ApplyTransformation (const Transform3D& t)
{
    if(!t.isIdentity())
        PSvectorTransform3D(t).Apply(pArray);
    return true;
}

void ParticleBunch::SortByCT ()
{
    //	pArray.sort();
    SortArray(pArray);
}

void ParticleBunch::Output (std::ostream& os) const
{
    //	std::copy(begin(),end(),ostream_iterator<PSvector>(os));
    int oldp=os.precision(10);
    ios_base::fmtflags oflg = os.setf(ios::scientific,ios::floatfield);
    for(PSvectorArray::const_iterator p = begin(); p!=end(); p++) {
        os<<std::setw(24)<<GetReferenceTime();
        os<<std::setw(24)<<GetReferenceMomentum();
        for(size_t k=0; k<6; k++)
            os<<std::setw(20)<<(*p)[k];
        os<<endl;
    }
    os.precision(oldp);
    os.flags(oflg);
}

void ParticleBunch::SetCentroid (const Particle& x0)
{
    PSvector x;
    GetCentroid(x);
    x-=x0;
    for(PSvectorArray::iterator p = begin(); p!=end(); p++)
        *p-=x;
}

}; // end namespace ParticleTracking

