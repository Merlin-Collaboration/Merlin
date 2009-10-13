// class RTMap
//
// A second-order TRANSPORT map for a PSvector
//
// class RTMap represents both the first-order R and second-order T TRANSPORT
// matrices. For efficiency, only non-zero terms are stored.

#ifndef RMap_h
#define RMap_h 1

#include "merlin_config.h"
#include <cassert>
#include <vector>
#include "BeamModel/PSTypes.h"
#include "TLAS/LinearAlgebra.h"
#include "NumericalUtils/utils.h"

// class RMap
// A linear phase space map. RMap represents a 6x6 matrixs (R matrix) which
// can map either a single PSvector or a SigmaMatrix (2nd-order moments), or
// a PSmoment object (combination of first- and second-order moments).
// Efficiency in terms of space and speed are achieved by using a sparse-matrix
// representation.

class RMap  {
public:

    // A single term in the matrix
    struct Rij {
        int i,j;
        double val;
        Rij(){}
        Rij(int k, int l, double v =0) : i(k),j(l),val(v) {}
        void Apply(const PSvector& orig,PSvector& res) const {
            res[i]+=val*orig[j];
        }
    };

    //	typedef std::vector<Rij> LinearTermArray;
    //	typedef LinearTermArray::iterator itor;
    //	typedef LinearTermArray::const_iterator const_itor;

    struct LinearTermArray {

        typedef Rij* iterator;
        typedef const Rij* const_iterator;

        Rij array[36];
        Rij* last;
        LinearTermArray() {last=array;}
        void push_back(const Rij& r) { *last = r; last++; }
        Rij& back() {return *(last-1);}
        const Rij& back() const {return *(last-1);}
        void reserve(size_t) {/* dummy */}
        Rij* begin() { return array; }
        Rij* end() { return last; }
        const Rij* begin() const { return array; }
        const Rij* end() const { return last; }
    };

    typedef LinearTermArray::iterator itor;
    typedef LinearTermArray::const_iterator const_itor;

public:

    RMap() : rterms() { /*rterms.reserve(12);*/ }
    explicit RMap(const RealMatrix&);

    // Mapping functions

    // Apply to a single PSvector
    PSvector& Apply(PSvector&) const;

    // Apply to a covariance (sigma) matrix
    // - MSVC6 needs this to be defined here (inlined) for
    // - the compiler to correctly compile it!
    template<class T, int N>
    TCovMtrx<T,N>& Apply(TCovMtrx<T,N>& S) const {
        TCovMtrx<T,N> result;
        for(const_itor t1 = rterms.begin(); t1!=rterms.end(); t1++) {
            for(const_itor t2 = rterms.begin(); t2!=rterms.end(); t2++) {
                map2nd(result,*t1,*t2,S);
            }
        }
        return S=result;
    }

    // Apply to PSmoments
    template<int N>
    TPSMoments<N>& Apply(TPSMoments<N>& S) const {
        TPSMoments<N> result;
        for(const_itor t1 = rterms.begin(); t1!=rterms.end(); t1++) {
            t1->Apply(S,result);
            for(const_itor t2 = rterms.begin(); t2!=rterms.end(); t2++) {
                map2nd(result,*t1,*t2,S);
            }
        }
        return S=result;
    }

    // Conversion to a matrix
    void ToMatrix(RealMatrix&, bool init=true) const;
    operator RealMatrix () const;

    // Output (matrix form)
    void MatrixForm(std::ostream&) const;

    // Array like accessors
    double operator()(int i, int j) const;
    double& operator()(int i, int j);

    // Term construction
    void AddTerm(int i, int j, double v) {
        rterms.push_back(Rij(i-1,j-1,v));
    }
    double& AddTerm(int i, int j) {
        rterms.push_back(Rij(i-1,j-1));
        return rterms.back().val;
    }

protected:

    void Apply(const PSvector&, PSvector&) const;

    // helper function used to map second order moments
    template<class T, int N>
    void map2nd(TCovMtrx<T,N>& res, const Rij& r1, const Rij& r2,
                const TCovMtrx<T,N>& orig) const {
        if(r1.i<N && r2.i<N && r1.j<N && r2.j<N) {
            double v = r1.val*r2.val*orig(r1.j,r2.j);
            res(r1.i,r2.i)+=(r1.i==r2.i) ? v : v/2.0;
        }
    }

    itor FindTerm(int,int);

    // non-zero terms
    LinearTermArray rterms;
};

// class R2Map
// simple POD R matrix representation for one degree of freedom
struct R2Map {
    R2Map() {}
    R2Map(double a11, double a12, double a21, double a22)
            :r11(a11),r12(a12),r21(a21),r22(a22) {}

    // array subscript operators (i={1,2})
    double* operator[](int n) { return &r11+2*n-3; }
    const double* operator[](int n) const { return &r11+2*n-3; }

    bool IsIdentity() const { return fequal(r11,1.0,1.0e-6) && fequal(r22,1.0,1.0e-6); }

    // matrix elements
    double r11,r12,r21,r22;
};

// Routines for applying R2Map objects

// Applies the R2Map to the specified 'plane' of a PSvector
inline void ApplyR2Map(const R2Map& R, PSvector& X, int plane)
{
    double *x0 = &X[0]+2*plane;
    const double x1 = R.r11*x0[0]+R.r12*x0[1];
    const double x2 = R.r21*x0[0]+R.r22*x0[1];
    x0[0]=x1;
    x0[1]=x2;
}

// Map (transform) plane terms: R.S.R'
inline void ApplyR2Map(const R2Map& R, TPSMoments<2>& S, int plane)
{
    ApplyR2Map(R,static_cast<PSvector&>(S),plane);

    const int p = 2*plane;

    double s11 = R.r11*R.r11*S(p,p);
    s11+=(R.r12*R.r11+R.r11*R.r12)*S(p+1,p);
    s11+=R.r12*R.r12*S(p+1,p+1);

    double s12 = R.r11*R.r21*S(p,p);
    s12+=(R.r12*R.r21+R.r11*R.r22)*S(p+1,p);
    s12+=R.r12*R.r22*S(p+1,p+1);

    double s22 = R.r21*R.r21*S(p,p);
    s22+=(R.r22*R.r21+R.r21*R.r22)*S(p+1,p);
    s22+=R.r22*R.r22*S(p+1,p+1);

    S(p,p)=s11;
    S(p+1,p)=s12;
    S(p+1,p+1)=s22;
}

// Map (transform) coupled terms Rx.Sxy.Ry'
inline void ApplyR2Map(const R2Map& Rx, TPSMoments<2>& S, const R2Map& Ry)
{
    double s13=0;
    s13+=Rx.r11*Ry.r11*S(0,2);
    s13+=Rx.r11*Ry.r12*S(0,3);
    s13+=Rx.r12*Ry.r11*S(1,2);
    s13+=Rx.r12*Ry.r12*S(1,3);

    double s14=0;
    s14+=Rx.r11*Ry.r21*S(0,2);
    s14+=Rx.r11*Ry.r22*S(0,3);
    s14+=Rx.r12*Ry.r21*S(1,2);
    s14+=Rx.r12*Ry.r22*S(1,3);

    double s23=0;
    s23+=Rx.r21*Ry.r11*S(0,2);
    s23+=Rx.r21*Ry.r12*S(0,3);
    s23+=Rx.r22*Ry.r11*S(1,2);
    s23+=Rx.r22*Ry.r12*S(1,3);

    double s24=0;
    s24+=Rx.r21*Ry.r21*S(0,2);
    s24+=Rx.r21*Ry.r22*S(0,3);
    s24+=Rx.r22*Ry.r21*S(1,2);
    s24+=Rx.r22*Ry.r22*S(1,3);

    S(0,2)=s13;
    S(0,3)=s14;
    S(1,2)=s23;
    S(1,3)=s24;
}


// global functions
void MakeIdentity(RMap& R, int ndf=3);

// utility routines for mapping arrays of objects

// Functor for applying a map as given
template<class M,class X>
struct map_applicator {
    const M& m;
    map_applicator(const M& anM) : m(anM) {}
    void operator()(X& x) const {
        m.Apply(x);
    }
};

// Functions for apply simple drift spaces.
// We provide these functions for effeciency: since drift
// spaces generally make up half the number of elements
// in a beamline, we provide special 'tuned' routines
// for speeding the drifts up.

// Apply a simple drift
struct ApplySimpleDrift {
    double z;

    // Apply to single vector
    explicit ApplySimpleDrift(double len):z(len) {}
    void Apply(PSvector& x) const {
        x.x()+=z*x.xp();
        x.y()+=z*x.yp();
    }

    // apply to 4D phase space moments
    void Apply(TPSMoments<2>& s) const {
        Apply(static_cast<PSvector&>(s));
        s(0,0)+=z*(2*s(0,1)+z*s(1,1));
        s(0,1)+=z*s(1,1);
        s(0,2)+=z*(s(0,3)+s(1,2)+z*s(1,3));
        s(0,3)+=z*s(1,3);
        s(1,2)+=z*s(1,3);
        s(2,2)+=z*(2*s(2,3)+z*s(3,3));
        s(2,3)+=z*s(3,3);
    }

    // apply to 6D phase space moments
    void Apply(TPSMoments<3>& s) const {
        Apply(static_cast<PSvector&>(s));
        s(0,0)+=z*(2*s(0,1)+z*s(1,1));
        s(0,1)+=z*s(1,1);
        s(0,2)+=z*(s(0,3)+s(1,2)+z*s(1,3));
        s(0,3)+=z*s(1,3);
        s(1,2)+=z*s(1,3);
        s(2,2)+=z*(2*s(2,3)+z*s(3,3));
        s(2,3)+=z*s(3,3);
        s(0,4)+=z*s(1,4);
        s(0,5)+=z*s(1,5);
        s(2,4)+=z*s(3,4);
        s(2,5)+=z*s(3,5);
    }

    template<class T>
    void operator()(T& s) const { Apply(s); }
};

// Apply a simple drift with (2nd-order) path length
struct ApplyDriftWithPathLength {
    double s;
    explicit ApplyDriftWithPathLength(double len):s(len) {}
    void Apply(PSvector& x) const {
        x.x()+=s*x.xp();
        x.y()+=s*x.yp();
        x.ct()-=s*(x.xp()*x.xp()+x.yp()*x.yp())/2.0;
    }
    void operator()(PSvector& x) const {
        Apply(x);
    }
};

// Functor for apply a map after re-scaling the momentum.
// Used mainly for tracking dipole magnets where the
// reference geometry (curvature) does not match the bend
// field / beam energy. typename X must supply double X::dp()
// function.
template<class M,class X>
struct map_applicator_dp {
    const M& m;
    double p_ratio;
    map_applicator_dp(const M& anM, double dp_r) : m(anM), p_ratio(dp_r) {}
    void operator()(X& x) const {
        double dp = x.dp();
        x.dp() = p_ratio*(1+dp)-1;
        m.Apply(x);
        x.dp() = dp;
    }
};

template<class C, class M>
void ApplyMap(const M& m, C& cont)
{
    for_each(cont.begin(),cont.end(),map_applicator<M,__TYPENAME__ C::value_type>(m));
}

template<class C, class M>
void ApplyMap(const M& m, C& cont, double p0, double p1)
{
    for_each(cont.begin(),cont.end(),map_applicator_dp<M,__TYPENAME__ C::value_type>(m,p0/p1));
}

#endif
