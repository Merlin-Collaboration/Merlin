/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

// Contains low-level BLAS-like routines
// for Matrix-Vector algebra.
// (Based on the FORTRAN BLAS library)

#ifndef _h_tblas
#define _h_tblas

#include <cassert>

namespace tblas1 {

// dot <- alpha + x'.y
template<class V1, class V2, class T>
T tdot(const V1& v1, const V2& v2, T alpha)
{
    for(register size_t i=0; i<v1.size(); ++i)
        alpha+=v1(i)*v2(i);
    return alpha;
}

template<class V1,class V2, class T>
void taxpy(T alpha, const V1& x, V2& y)
{
    assert(&x!=&y);
    if(alpha==0)
        return;
    if(alpha==1) {
        y+=x;
    }
    else {
        for(register size_t i=0; i<x.size(); ++i)
            y(i)+=alpha*x(i);
    }
}
};


namespace tblas2 {

using namespace tblas1;

// y <- alpha*A.x + beta*y
template<class Ta, class M, class Vx, class Tb, class Vy>
void tgemv(bool t, const Ta& alpha, const M& A, const Vx& x, const Tb& beta, Vy& y)
{
    assert(&x!=&y); // that's a no-no!
    if(beta!=1)
        for(register size_t i=0; i<A.nrows(); i++)
            y[i]*=beta;

    if(alpha!=0)
        if(t) { // use transpose of A
            for(register size_t i=0; i<A.ncols(); ++i)
                for(register size_t j=0; j<A.nrows(); ++j)
                    y[i]+=alpha*A(j,i)*x[j];
        }
        else { // use normal A
            for(register size_t i=0; i<A.nrows(); ++i)
                for(register size_t j=0; j<A.ncols(); ++j)
                    y[i]+=alpha*A(i,j)*x[j];
        }
}

// add Symm, Upper and Lower diagonal forms later

// A <- alpha*x.y' + A (outer product)
template<class T, class Vx, class Vy, class M>
void tger(const T& alpha, const Vx& x, const Vy& y, M& A)
{
    if(alpha==0)
        return;

    for(register size_t i=0; i<A.nrows(); ++i)
        for(register size_t j=0; j<A.ncols(); ++j)
            A(i,j)+=alpha*x(i)*y(j);
}
}

namespace tblas3 {

using namespace tblas2;

// C <- alpha*A.B + beta*C
template<class Ta, class Ma, class Mb, class Tb, class Mc>
void tgemm(bool tpa, bool tpb, const Ta& alpha, const Ma& A, const Mb& B,
           const Tb& beta, Mc& C)
{
    assert(&C!=&A && &C!=&B); // that's a no-no

    if(beta==0)
        C=0; // assignment faster than multiplication?
    else if(beta!=1)
        C*=beta;

    if(alpha==0)
        return;

    // calculate the range of the k subscript
    const size_t klim = tpa ? A.nrows() : A.ncols();

    for(register size_t i=0; i<C.nrows(); ++i)
        for(register size_t j=0; j<C.ncols(); ++j)
            for(register size_t k=0; k<klim; ++k) {
                // 4 cases to deal with here
                if(tpa) {
                    if(tpb)
                        C(i,j)+=alpha*A(k,i)*B(j,k); // A'.B'
                    else
                        C(i,j)+=alpha*A(k,i)*B(k,j); // A'.B
                }
                else {
                    if(tpb)
                        C(i,j)+=alpha*A(i,k)*B(j,k); // A.B'
                    else
                        C(i,j)+=alpha*A(i,k)*B(k,j); // A.B
                }
            }
}

// Matrix rotation C <- R.M.R'
template<class Tr, class Ta, class Tb>
void tgemr(const Tr& R, const Ta& M, Tb& C)
{
    const size_t n = R.nrows();
    for(register size_t i=0; i<n; i++)
        for(register size_t j=0; j<n; j++) {
            C(i,j)=0;
            for(register size_t k=0; k<n; k++)
                for(register size_t l=0; l<n; l++)
                    C(i,j)+=R(i,k)*R(j,l)*M(k,l);
        }
}
}

#endif
