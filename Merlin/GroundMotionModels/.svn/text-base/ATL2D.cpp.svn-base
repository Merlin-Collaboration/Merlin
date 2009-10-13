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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <algorithm>
#include <iomanip>

// RandomNG
#include "Random/RandomNG.h"
// ATL2D
#include "GroundMotionModels/ATL2D.h"

#undef ATL_XY

namespace {

using namespace std;

inline void ResetSupport(AcceleratorSupport* s)
{
    s->Reset();
}

struct ApplyATL {

    ApplyATL(RealVector& gmy, double vibv, RandGenerator* rg, ATL2D::ATLMode theMode)
            : n(0),yy(gmy),vv(vibv),rng(rg),atlMode(theMode) {}

    void operator()(AcceleratorSupport* s)
    {
        if(atlMode == ATL2D::increment)
            s->IncrementOffset(0,yy(n),0);
        else
        {
            // add random 'noise'
            double yv = vv!=0 ? rng->normal(0,vv) : 0.0;
            s->SetOffset(0,yy(n)+yv,0);
        }

        n++;
    }

    // data members
    int n;
    RealVector yy;
    double vv;
    RandGenerator* rng;
    ATL2D::ATLMode atlMode;
};

struct DumpOffset {

    DumpOffset(ostream& anOS, const double t): os(anOS),time(t) {}

    void operator()(const AcceleratorSupport* s)
    {
        using std::setw;
        os<<setw(14)<<time;
        double arcs =s->GetArcPosition();
        os<<setw(14)<<arcs;
        Point2D locn=s->GetLocation();
        os<<setw(14)<<locn.x;
        os<<setw(14)<<locn.y;
        Vector3D offset=s->GetOffset();
        //			os<<setw(14)<<offset.x;
        os<<setw(14)<<offset.y;
        //			os<<setw(14)<<offset.z;
        os<<endl;
    }

    ostream& os;
    double time;
};

}; // end of anonymous namespace


ATL2D::ATL2D (double anA, const AcceleratorSupportList& supports, const Point2D refPoint, ifstream* evecTFile, ifstream* evalFile)
        : t(0),A(anA),seed(0),vv(0),theSupports(supports),rg(new RandGenerator()),atlMode(absolute)
{
    const int n = theSupports.size();

    evecsT.redim(n,n);
    evals.redim(n);
    int row,col;

    if(evecTFile && evalFile)
    {
        double element;
        for(row=0;   row<n; row++)
            for(col=row; col<n; col++)
            {
                (*evecTFile)>>element;
                evecsT(row,col) = evecsT(col,row) = element;
            }

        for(row=0; row<n; row++)
            (*evalFile)>>evals(row);
    }
    else
    {
        for(row=0;   row<n; row++)
            for(col=row; col<n; col++)
                evecsT(row,col) = evecsT(col,row) = (Distance(col,refPoint) + Distance(row,refPoint) - Distance(row,col))/2.0;

        EigenSystemSymmetricMatrix(evecsT, evals);
    }

    rg->init(seed);
}

ATL2D::~ATL2D ()
{
    delete rg;
}

void ATL2D::Reset ()
{
    for_each(theSupports.begin(),theSupports.end(),ResetSupport);
    t=0;
}

double ATL2D::DoStep (double dt)
{
    RealVector yy(theSupports.size());

    double at = (atlMode==increment) ? A * dt : A * t;

    for(int n=0; n<theSupports.size(); n++)
        yy(n) = rg->normal(0,at*evals[n]);

    RealVector dy = evecsT * yy;

    for_each(theSupports.begin(),theSupports.end(),ApplyATL(dy,vv,rg,atlMode));

    return t+=dt;
}

void ATL2D::RecordOffsets (std::ostream& os) const
{
    ios_base::fmtflags oldFlags = os.flags();
    int prec = os.precision();

    os.setf(ios_base::scientific,ios_base::floatfield);
    os.precision(4);

    for_each(theSupports.begin(),theSupports.end(),DumpOffset(os,t));

    os.flags(oldFlags);
    os.precision(prec);
}

double ATL2D::GetTime () const
{
    return t;
}

void ATL2D::SetRandomSeed (unsigned int nseed)
{
    rg->reset(nseed);
}

unsigned int ATL2D::GetRandomSeed () const
{
    return rg->getSeed();
}

void ATL2D::ResetRandomSeed ()
{
    rg->reset();
}

bool ATL2D::SetATLMode(const ATLMode mode)
{
    if( (atlMode==increment) && (vv>0) )
        return false;

    atlMode = mode;
    return true;
}

bool ATL2D::SetVibration(const double vrms)
{
    if(atlMode==increment)
        return false;

    vv = vrms*vrms;
    return true;
}

void ATL2D::RecordEigenSystem(ofstream* evecTFile, ofstream* evalFile)
{
    const int n = evals.size();
    int row,col;

    for(row=0;   row<n; row++)
    {
        for(col=row; col<n; col++)
            (*evecTFile)<<std::setw(14)<<evecsT(row,col);
        (*evecTFile)<<endl;
    }

    for(row=0; row<n; row++)
        (*evalFile)<<std::setw(14)<<evals(row)<<endl;
}

double ATL2D::Distance(const int n1, const int n2)
{
    const Point2D x1 = theSupports[n1]->GetLocation();
    const Point2D x2 = theSupports[n2]->GetLocation();

    const Point2D dx = x1 - x2;

    return sqrt( dx.x * dx.x + dx.y * dx.y );
}

double ATL2D::Distance(const int n1, const Point2D x2)
{
    const Point2D x1 = theSupports[n1]->GetLocation();

    const Point2D dx = x1 - x2;

    return sqrt( dx.x * dx.x + dx.y * dx.y );
}

