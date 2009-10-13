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
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_TLASimp
#define _h_TLASimp

// Templated functions and classes for performing linear algebra. The following
// perform linear algebra functions based on TMAT::Matrix and TMAT::Vector classes.

#include <vector>
#include "TLAS/TLAS.h"

namespace TLAS {
// Class Declarations

// LUMatrix
// Represents a LU-decomposition of a square matrix. An LUMatrix object
// can be used to "solve" a system of linear equations with many RHS vectors.

template<class T>
class LUMatrix {
public:

    // Construction from an arbitrary (square) matrix.
    // Throws DimensionError() if the matrix is not square, or
    // SingularMatrix() if the matrix is singular.
    explicit LUMatrix(const Matrix<T>& M) throw (DimensionError, SingularMatrix)
            : lud(M),indecies(),d(1)
    {
        ludcmp(lud,indecies,d);
    }

    // Solve the RHS vector using this LU-decomposition.
    // Note that rhs is over-written with the result.
    // Returns rhs.
    Vector<T>& operator()(Vector<T>& rhs) const {
        return lubksb(lud,indecies,rhs);
    }
    SubVector<T>& operator()(SubVector<T>& rhs) const {
        return lubksb(lud,indecies,rhs);
    }

    // These functions are provided for acting on const. rhs vectors.
    Vector<T> operator()(const Vector<T>& rhs) const {
        Vector<T> v(rhs);
        return operator()(v);
    }
    Vector<T> operator()(const SubVector<T>& rhs) const {
        Vector<T> v(rhs);
        return operator()(v);
    }
    Vector<T> operator()(const ConstSubVector<T>& rhs) const {
        Vector<T> v(rhs);
        return operator()(v);
    }

    // Return the determinant of the original matrix
    T det() const {
        double dt=d;
        for(int i=0;i<lud.ncols();i++)
            dt*=lud(i,i);
        return dt;
    }

private:
    Matrix<T> lud;
    std::vector<int> indecies;
    T d;
};

// SVDMatrix
// Represents a Singular-Value Decomposition of a matrix. As with LUMatrix, an
// SVDMatrix object can be used to solve a system of linear equations with many
// RHS vectors.

template<class T>
class SVDMatrix {
public:

    // Construction from an arbitrary matrix. If M.nrows()<M.ncols(), then M is first
    // made square by the addition of M.ncols()-M.nrows() zero rows. threshold specifies
    // the relative (to the largest singular value) threshold value below which the
    // singular values are set to zero.
    explicit SVDMatrix(const Matrix<T>& M, T threshold = T(1e-06)) throw (ConvergenceFailure,SingularValuesAllZero);
    SVDMatrix(const Matrix<T>& M, const Vector<T>& wts, T threshold = T(1e-06)) throw (ConvergenceFailure,SingularValuesAllZero);

    // Solve the RHS vector using this SVD.
    Vector<T> operator()(const Vector<T>& rhs) const {
        Vector<T> x(w.size());
        Vector<T> y(rhs);
        y*=wts;
        return svbksb(u,w,v,y,x);
    }
    Vector<T> operator()(const SubVector<T>& rhs) const {
        Vector<T> x(w.size());
        Vector<T> y(rhs);
        y*=wts;
        return svbksb(u,w,v,y,x);
    }
    Vector<T> operator()(const ConstSubVector<T>& rhs) const {
        Vector<T> x(w.size());
        Vector<T> y(rhs);
        y*=wts;
        return svbksb(u,w,v,y,x);
    }
    const std::vector<bool>& wflags() const { return wflgs; }

    // Decomposed matrices
    const Matrix<T>& U() const { return u; }
    const Matrix<T>& V() const { return v; }
    const Vector<T>& W() const { return w; }

private:

    void Init(const Matrix<T>& M, T threshold);

    Matrix<T> u,v;
    Vector<T> w;
    Vector<T> wts;
    std::vector<bool> wflgs;
};

template<class T>
SVDMatrix<T>::SVDMatrix(const Matrix<T>& M, T threshold) throw (ConvergenceFailure, SingularValuesAllZero)
        : wts(M.nrows())
{
    wts=T(1);
    Init(M,threshold);
}

template<class T>
SVDMatrix<T>::SVDMatrix(const Matrix<T>& M, const Vector<T>& wts1, T threshold)
throw (ConvergenceFailure, SingularValuesAllZero)
        : wts(M.nrows())
{
    wts=wts1;
    Init(M,threshold);
}

template<class T>
void SVDMatrix<T>::Init(const Matrix<T>& M, T threshold)
{
    if(M.nrows()<M.ncols()) {
        u.redim(M.ncols(),M.ncols());
        u(Range(0,M.nrows()-1),Range(0,M.ncols()-1))=M;
    }
    else
        u.copy(M);

    // adjust for weights: multiply each column vector by wts
    for(Subscript ir=0; ir<u.ncols(); ir++)
        u.column(ir) *= wts;

    w.redim(u.ncols());
    wflgs = std::vector<bool>(u.ncols(),true);

    v.redim(u.ncols(),u.ncols());
    svdcmp(u,w,v);

    T wmin = threshold!=T(0) ? threshold*(*max_element(w.begin(),w.end())) : threshold;
    int zerocount=0;
    for(size_t i=0; i<w.size(); i++)
        if(w[i]<=wmin) {
            w[i]=T(0);
            wflgs[i]=false;
            zerocount++;
        }
    if(zerocount==static_cast<int>(w.size())) {
        throw SingularValuesAllZero();
    }
}

// Low level routines (modified from Numerical Recipes in C)

template<class T>
void ludcmp(Matrix<T>& a, std::vector<int>& indecies, T& d)
{
    static const T tiny = numeric_limits<T>::epsilon();
    int imax,i,j,k;
    const int n = a.ncols();

    if(n!=a.nrows())
        throw DimensionError();

    std::vector<T> vv(n);
    std::vector<int> indx(n);

    d=1.0;
    for (i=0;i<n;i++) {
        double big=0.0;
        for (j=0;j<n;j++)
            big = _MAX(fabs(a[i][j]),big);
        if (big == 0.0)
            throw SingularMatrix();

        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            double sum=a[i][j];
            for (k=0;k<i;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        double big=0.0;
        for (i=j;i<n;i++) {
            double sum=a[i][j];
            for (k=0;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            double dum;
            if ((dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0;k<n;k++) {
                double dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            d*=-1;
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0)
            a[j][j] = tiny;
        if (j<n) {
            double dum=1.0/(a[j][j]);
            for (i=j+1;i<n;i++)
                a[i][j] *= dum;
        }
    }

    // copy the indecies
    indecies.swap(indx);
}

// lu back substitution

template<class T, class V>
V& lubksb(const Matrix<T>& a, const std::vector<int>& indx, V& b)
{
    int i,ii=-1,ip,j;
    const int n = a.ncols();
    T sum;

    for(i=0;i<n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if(ii!=-1) {
            for (j=ii;j<=i-1;j++)
                sum -= a[i][j]*b[j];
        }
        else if(sum!=0)
            ii=i;
        b[i]=sum;
    }
    for (i=n-1;i>=0;i--) {
        sum=b[i];
        for (j=i+1;j<n;j++)
            sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
    return b;
}

template<class T>
Matrix<T>& InvertMatrix(Matrix<T>& m)
{
    if(m.nrows()!=m.ncols())
        throw NonSquareMatrix();

    LUMatrix<T> lu(m);
    for(int i = 0; i<m.ncols(); i++) {
        SubVector<T>& col = m.column(i);
        col=T(0);
        col[i]=T(1);
        lu(col);
    }
    return m;
}

template<class T> inline Matrix<T> Invert(const Matrix<T>& m)
{
    Matrix<T> local(m);
    return InvertMatrix(local);
}

template<class T> inline T Det(const Matrix<T>& m)
{
    return LUMatrix<T>(m).det();
}

// SVD routines
	#include "TLAS/svdcmp.h"

}; // end namespace TLAS

#endif // _h_TLASimp
