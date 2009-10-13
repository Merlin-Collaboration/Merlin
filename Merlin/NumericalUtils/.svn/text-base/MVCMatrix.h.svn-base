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

#ifndef MCVMatrix_h
#define MCVMatrix_h 1

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

#define MVCMTRX_REDUCE_MEMORY 1

using std::ostream;
using std::setprecision;
using std::setw;
using std::right;
using std::fixed;
using std::scientific;
using std::endl;

inline int ComputeIndex(int i, int j)
{
    if(i<j) std::swap(i,j);
    return i*(i+1)/2+j;
}

inline int ComputeIndex(int i)
{
    return i*(i+3)/2;
}


//	An NxN matrix representing the full covariance matrix of
//	an N-variate data sample (including first-order moments,
//	or means). Variables are indexed from 0-(n-1).

template <class T, int N>
class MVCMatrix
{
public:
    //	Default constructor initialises all values to zero.
    MVCMatrix ();


    //	Return the mean of the n-th variable.
    T mean (int n) const;

    //	Return the standard diviation of the n-th variable.
    T std (int n) const;

    //	Return the variance of the n-th variable.
    T variance (int n) const;

    //	Return the correlation coefficient of the n-th and m-th
    //	variables =  s_mn/sqrt(s_mm*s_nn)
    T r (int n, int m) const;

    //	Matrix indexing. Returns <xi*xj>. Indexing runs from 0
    //	to N-1.
    T& operator () (int i, int j);

    const T& operator () (int i, int j) const;

    //	Returns the mean value of the n-th parameter.
    T& operator [] (int n);

    const T& operator [] (int n) const;

    //	Sets all elements to zero.
    void zero ();

    //	Prints the statistical data in a formatted table. If
    //	normalise==true, then diagonal terms are printed as
    //	correlation coefficients rij.
    void printFormatted (ostream& os, bool normalised = true) const;

    // logical comparisons
    bool operator==(const MVCMatrix<T,N>& rhs) const;
    bool operator!=(const MVCMatrix<T,N>& rhs) const {
        return !operator==(rhs);
    }

protected:
private:
    // Data Members for Class Attributes

    //	An array of the mean values of the distribution.
    T m1[N];

    //	An array of the second-order moments (variances) of the
    //	distribution.
    T m2[N*(1+N)/2];

private:
};



// Parameterized Class MVCMatrix


template <class T, int N>
inline T MVCMatrix<T,N>::mean (int n) const
{
    return m1[n];
}

template <class T, int N>
inline T MVCMatrix<T,N>::std (int n) const
{
#ifndef MVCMTRX_REDUCE_MEMORY
    return sqrt(m2[n][n]);
#else
    return sqrt(m2[ComputeIndex(n)]);
#endif
}

template <class T, int N>
inline T MVCMatrix<T,N>::variance (int n) const
{
#ifndef MVCMTRX_REDUCE_MEMORY
    return m2[n][n];
#else
    return m2[ComputeIndex(n)];
#endif
}

template <class T, int N>
inline T MVCMatrix<T,N>::r (int n, int m) const
{
    return n==m ? 1 :
#ifndef MVCMTRX_REDUCE_MEMORY
           m2[n][m]/sqrt(variance(n)*variance(m));
#else
           m2[ComputeIndex(n,m)]/sqrt(variance(n)*variance(m));
#endif
}

template <class T, int N>
inline T& MVCMatrix<T,N>::operator () (int i, int j)
{
#ifndef MVCMTRX_REDUCE_MEMORY
    return m2[i][j];
#else
    return m2[ComputeIndex(i,j)];
#endif
}

template <class T, int N>
inline const T& MVCMatrix<T,N>::operator () (int i, int j) const
{
#ifndef MVCMTRX_REDUCE_MEMORY
    return m2[i][j];
#else
    return m2[ComputeIndex(i,j)];
#endif
}

template <class T, int N>
inline T& MVCMatrix<T,N>::operator [] (int n)
{
    return m1[n];
}

template <class T, int N>
inline const T& MVCMatrix<T,N>::operator [] (int n) const
{
    return m1[n];
}

template <class T, int N>
inline void MVCMatrix<T,N>::zero ()
{
#ifndef MVCMTRX_REDUCE_MEMORY
    std::fill(m1,reinterpret_cast<T*>(m2)+N*N,T(0));
#else
    std::fill(m1,m2+N*(1+N)/2,T(0));
#endif
}





// Parameterized Class MVCMatrix



template <class T, int N>
MVCMatrix<T,N>::MVCMatrix ()
{
    zero();
}



template <class T, int N>
void MVCMatrix<T,N>::printFormatted (ostream& os, bool normalised) const
{
    using std::setprecision;
    using std::setw;
    for(int i=0; i<N; i++) {
        os<<setw(12)<<right<<scientific<<setprecision(3)<<mean(i);
        os<<setw(12)<<right<<scientific<<setprecision(3)<<std(i);
        for(int j=1; j<N; j++) {
            if(j<=i) {
                if(normalised)
                    os<<setw(8)<<" ";
                else
                    os<<setw(12)<<" ";
            }
            else {
                os<<right<<setprecision(3);
                if(normalised)
                    os<<setw(8)<<fixed<<r(i,j);
                else
                    os<<setw(12)<<scientific<<operator()(i,j);
            }
        }
        os<<endl;
    }
}

template<class T, int N>
bool MVCMatrix<T,N>::operator==(const MVCMatrix<T,N>& rhs) const
{
    register int i;
    for(i=0;i<N;i++)
        if(m1[i]!=rhs.m1[i])
            return false;

    for(i=0;i<N*(1+N)/2;i++)
        if(m2[i]!=rhs.m2[i])
            return false;
    return true;
}

#endif
