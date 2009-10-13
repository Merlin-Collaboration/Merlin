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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef PSvector_h
#define PSvector_h 1

#include "merlin_config.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include "TLAS/LinearAlgebra.h"

// Phase space coordinate indecies
typedef int PScoord;
#define ps_X  0
#define ps_XP 1
#define ps_Y  2
#define ps_YP 3
#define ps_CT 4
#define ps_DP 5

class PSvector {
public:

    PSvector ()
    {}

    explicit PSvector (double x) {
        std::fill(v,v+6,x);
    }

    //	Component accessors.
    double x ()  const{ return v[0];}
    double y ()  const{ return v[2];}
    double ct () const{ return v[4];}
    double xp () const{ return v[1];}
    double yp () const{ return v[3];}
    double dp () const{ return v[5];}

    //	Array access.
    double operator [] (PScoord coord) const{ return v[coord];}

    //	Component mutators.
    double& x () { return v[0];}
    double& y () { return v[2];}
    double& ct (){ return v[4];}
    double& xp (){ return v[1];}
    double& yp (){ return v[3];}
    double& dp (){ return v[5];}

    //	Array access.
    double& operator [] (PScoord coord) {return v[coord];}

    //	Conversion to a RealVector.
    operator RealVector () const {return RealVector(v,6);}

    bool operator == (const PSvector& psv) const{
        return memcmp(v,psv.v,6*sizeof(double))==0;
    }

    bool operator != (const PSvector& psv) const{
        return memcmp(v,psv.v,6*sizeof(double))!=0;
    }

    //	Sets the vector to zero.
    void zero () {
        std::fill(v,v+6,0.0);
    }

    //	Arithmetic assignment
    PSvector& operator += (const PSvector& p) {
        double *q=v;
        const double *r=p.v;
        while(q!=(v+6)) *(q++) += *(r++);
        return *this;
    }

    PSvector& operator -= (const PSvector& p) {
        double *q=v;
        const double *r=p.v;
        while(q!=(v+6)) *(q++) -= *(r++);
        return *this;
    }

    PSvector& operator *= (double x) {
        for(double *q = v; q!=v+6;q++)
            (*q)*=x;
        return *this;
    }

    PSvector& operator /= (double x) {
        for(double *q = v; q!=v+6;q++)
            (*q)/=x;
        return *this;
    }

    // binary arithmetic operators
    PSvector operator+(const PSvector& rhs) const {
        PSvector rv(*this);
        return rv+=rhs;
    }

    PSvector operator-(const PSvector& rhs) const {
        PSvector rv(*this);
        return rv-=rhs;
    }

    //	Reads the vector as six floating point numbers,
    //	separated by spaces and terminated by a newline.
    friend std::ostream& operator<<(std::ostream& os, const PSvector& v);

    //	Outputs the vector in row form, space delimeted with a
    //	terminating newline.
    friend std::istream& operator>>(std::istream& is, PSvector& v);

private:
    double v[6];
};

//	A linear array of PSvector objects.
typedef std::vector<PSvector> PSvectorArray;

#endif
