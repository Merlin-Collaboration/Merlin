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
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#include "NumericalUtils/Interpolation.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <cassert>

namespace {

using namespace std;

// Method used for arbitrarily space data
//
class ArbSpacedData : public Interpolation::Method {
public:

    struct Data {
        double x,y;
        Data(double x1=0, double y1=0) : x(x1),y(y1) {}
    };

    ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals);
    virtual double ValueAt(double) const throw(Interpolation::BadRange);

private:
    vector<Data> itsData;
};

// Method used for equally space data
//
class EqualSpacedData : public Interpolation::Method {
public:

    EqualSpacedData(const vector<double> yv, double xm, double delta)
            : yvals(yv),xmin(xm),xmax(xm+(yv.size()-1)*delta),dx(delta)
    {
        assert(dx>0);
    }

    double ValueAt(double x) const throw(Interpolation::BadRange);

private:
    vector<double> yvals;
    double xmin;
    double xmax;
    double dx;
};

// Function  used by sort() for ArbSpacedData::Data
//
inline bool sort_x(const ArbSpacedData::Data& d1, const ArbSpacedData::Data& d2)
{
    return d1.x<d2.x;
}

// Class ArbSpacedData implementation
//
ArbSpacedData::ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals)
        : itsData()
{
    assert(xvals.size()==yvals.size());
    itsData.reserve(xvals.size());
    for(size_t i=0; i<xvals.size(); i++)
        itsData.push_back(Data(xvals[i],yvals[i]));
    sort(itsData.begin(),itsData.end(),sort_x);
}

double ArbSpacedData::ValueAt(double x) const throw (Interpolation::BadRange)
{
    if(x<itsData.front().x || x>itsData.back().x)
        throw Interpolation::BadRange(x,FloatRange(itsData.front().x,itsData.back().x));

    // locate segment by binary search
    size_t ju=itsData.size(),jl=0,jm;

    while((ju-jl)>1) {
        jm = (ju+jl)>>1;
        if(x>itsData[jm].x)
            jl=jm;
        else
            ju=jm;
    }

    // use linear interpolation to return value
    double m = (itsData[jl+1].y-itsData[jl].y)/(itsData[jl+1].x-itsData[jl].x);
    return itsData[jl].y+(x-itsData[jl].x)*m;
}

// Class EqualSpacedData implementation
//

double EqualSpacedData::ValueAt(double x) const throw(Interpolation::BadRange)
{
    // note that we use extrapolation here if x is out of range
    size_t n;
    if(x<xmin)
        n=0;
    else if(x>xmax)
        n=yvals.size()-2;
    else
        n = (x-xmin)/dx;

    double m = (yvals[n+1]-yvals[n])/dx;
    double x0 = xmin+dx*n;
    return yvals[n]+m*(x-x0);
}

};

// exception
Interpolation::BadRange::BadRange(double x, const FloatRange& r)
        : MerlinException(""),value(x),valid_range(r)
{
    ostringstream buf;
    buf<<value<<" not in interpolation range ("<<r.lower<<","<<r.upper<<')';
    SetMsg(buf.str());
}

Interpolation::Interpolation(const vector<double>& yvals, double xmin, double dx)
        : itsMethod(new EqualSpacedData(yvals,xmin,dx))
{}

Interpolation::Interpolation(const std::vector<double>& xv, const vector<double>& yv)
        : itsMethod(new ArbSpacedData(xv,yv))
{}

Interpolation::~Interpolation()
{
    if(itsMethod)
        delete itsMethod;
}
