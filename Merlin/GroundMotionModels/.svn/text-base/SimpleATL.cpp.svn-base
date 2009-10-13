/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iomanip>
// RandomNG
#include "Random/RandomNG.h"
// SimpleATL
#include "GroundMotionModels/SimpleATL.h"

#undef ATL_XY

namespace {

using namespace std;

inline void ResetSupport(AcceleratorSupport* s)
{
    s->Reset();
}

// function used to sort supports in acending arc position
inline bool AsZ(const AcceleratorSupport* s1, const AcceleratorSupport* s2)
{
    return s1->GetArcPosition() < s2->GetArcPosition();
}

struct ApplyATL {

    ApplyATL(double A, double dt, vector<double>& gmy, double vibv,RandGenerator* rg)
            : AT(A*dt),y(0),z(0),yy(gmy.begin()),vv(vibv),rng(rg) {}

    void operator()(AcceleratorSupport* s)
    {
        // Perform the random walk
        double v = AT*(s->GetArcPosition()-z);
        y += rng->normal(0,v);
        *yy += y;

        // add random 'noise'
        double yv = vv!=0 ? rng->normal(0,vv) : 0.0;

        s->SetOffset(0,*yy+yv,0);
        z=s->GetArcPosition();

        yy++;
    }

    // data members
    double AT;
    double y;
    double z;
    vector<double>::iterator yy;
    double vv;
    RandGenerator* rng;
};

struct DumpOffset {

    DumpOffset(ostream& anOS):os(anOS) {}

    void operator()(const AcceleratorSupport* s)
    {
        using std::setw;
        Vector3D offset=s->GetOffset();
        os<<setw(12)<<offset.x;
        os<<setw(12)<<offset.y;
        os<<setw(12)<<offset.z;
        os<<endl;
    }

    ostream& os;
};

}; // end of annonymous namespace

SimpleATL::SimpleATL (double anA, const AcceleratorSupportList& supports, double vrms)
        : t(0),A(anA),seed(0),vv(vrms*vrms),atlgm(supports.size(),0.0),theSupports(supports),
          rg(new RandGenerator())
{
    sort(theSupports.begin(),theSupports.end(),AsZ);
}

SimpleATL::~SimpleATL ()
{
    delete rg;
}

void SimpleATL::Reset ()
{
    for_each(theSupports.begin(),theSupports.end(),ResetSupport);
    fill(atlgm.begin(),atlgm.end(),0.0);
    t=0;
}

double SimpleATL::DoStep (double dt)
{
    for_each(theSupports.begin(),theSupports.end(),ApplyATL(A,dt,atlgm,vv,rg));
    return t+=dt;
}

void SimpleATL::RecordOffsets (std::ostream& os) const
{
    ios_base::fmtflags oldFlags = os.flags();
    int prec = os.precision();

    os.setf(ios_base::scientific,ios_base::floatfield);
    os.precision(4);

    for_each(theSupports.begin(),theSupports.end(),DumpOffset(os));

    os.flags(oldFlags);
    os.precision(prec);
}

double SimpleATL::GetTime () const
{
    return t;
}

void SimpleATL::SetRandomSeed (unsigned int nseed)
{
    rg->reset(nseed);
}

unsigned int SimpleATL::GetRandomSeed () const
{
    return rg->getSeed();
}

void SimpleATL::ResetRandomSeed ()
{
    rg->reset();
}

