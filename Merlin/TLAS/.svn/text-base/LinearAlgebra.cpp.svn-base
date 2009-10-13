#include <fstream>
#include "TLAS/LinearAlgebra.h"
#include <algorithm> // for std::swap
#include <cmath>

#ifndef DEBUG
#include "NumericalUtils/MatrixPrinter.h"
#endif

namespace {

#ifndef DEBUG
ofstream debug_os("tlas_debug.dat");
#endif

// MSVC++ bug! Compiler should make correct resolution for abs()!!
inline double ABS(double x) { return fabs(x); }
inline double ABS(const Complex& z) { return abs(z); }

#define _SWAP(a,b) std::swap((a),(b))
#define _SIGN(a,b) ((b>0) ? fabs(a) : -fabs(a))

template<class T> double Inverse(Matrix<T>& a)
{

    if(a.nrows() != a.ncols())
        throw TLAS::NonSquareMatrix();

    int n = a.nrows();

    Vector<int> indxc(n);
    Vector<int> indxr(n);
    Vector<int> ipiv(n);

    int icol, irow,i,j,k;
    double minpiv=1.0e9;
    double big;
    T dum, pivinv;

    for(j=0; j<n; j++)
        ipiv(j)=0;

    for(i=0; i<n; i++)
    {
        big = 0.0;
        for(j=0; j<n; j++)
            if(ipiv(j) != 1)
                for(k=0; k<n; k++)
                {
                    if(ipiv(k)==0)
                    {
                        if(ABS(a(j,k))>=big)
                        {
                            big = ABS(a(j,k));
                            irow=j;
                            icol=k;
                        }
                    }
                    else if (ipiv[k]>1)
                        throw TLAS::SingularMatrix();
                }
        ++ipiv(icol);

        if(irow!=icol)
            for(int l=0; l<n; l++) _SWAP(a(irow,l),a(icol,l));

        indxr(i)=irow;
        indxc(i)=icol;
        if(ABS(a(icol,icol))<minpiv) minpiv=ABS(a(icol,icol));
        if(minpiv==0.0)
        {
            cout<<"Singular"<<endl;
            throw TLAS::SingularMatrix();
        }

        pivinv = 1.0 / a(icol,icol);
        a(icol,icol) = 1.0;
        for(int l=0; l<n; l++) a(icol,l) *= pivinv;

        for(int m=0; m<n; m++)
            if(m!=icol)
            {
                dum = a(m,icol);
                a(m,icol) = 0.0;
                for(int l=0; l<n; l++) a(m,l) -= a(icol,l)*dum;
            }
    }

    for(int l=n-1; l>=0; l--)
    {
        if(indxr(l)!=indxc(l))
            for(int k=0; k<n; k++) _SWAP(a(k,indxr(l)),a(k,indxc(l)));
    }
    return minpiv;
}
}; // end anonymous namespace

namespace TLAS {

double Invert(RealMatrix& t)
{
	return Inverse(t);
}

void EigenSystem(RealMatrix& t, ComplexVector& eigenvalues, ComplexMatrix& eigenvectors)
{
    ComplexMatrix m(t);
    Matrix<int> s(6,6,0);

    s(0,1) = s(2,3) = s(4,5) = 1;
    s(1,0) = s(3,2) = s(5,4) = -1;

    eigenvalues.redim(3);
    eigenvectors.redim(3,6);

    for(int pln=0; pln<3; pln++)
    {
        //Guess an eigenvalue
        Complex cos_mu((t(2*pln,2*pln)+t(2*pln+1,2*pln+1))/2.0,0);
        Complex isin_mu = sqrt(cos_mu*cos_mu - Complex(1,0));
        Complex lambda(cos_mu + isin_mu);

        //Guess an eigenvector
        ComplexVector b(0.0,6);
        b(2*pln) = Complex(1/sqrt(2.0));
        b(2*pln+1) = b(2*pln);

        //Iterate!
        ComplexMatrix mp(m);

        double prox = 1.0;
        int iter = 0;
        Subscript n=0;
        Subscript row=0;

        while(prox>1.0e-16 && iter<100)
        {
            for(n=0; n<6; n++)
                mp(n,n) = m(n,n) - lambda;

            ComplexMatrix minv(mp);
            prox = Inverse(minv);

            ComplexVector y(6);
            for(row=0; row<6; row++)
            {
                Complex sum=0.0;
                for(Subscript col=0; col<6; col++)
                    sum += minv(row,col) * b(col);
                y(row) = sum;
            }

            Complex ynorm = 0.0;
            Complex ydotb = 0.0;
            for(n=0; n<6; n++)
            {
                ynorm += y(n) * y(n);
                ydotb += y(n) * b(n);
            }

            b = y / sqrt(ynorm);
            lambda += 1.0/ydotb;
            iter++;
        }

        //Now normalise to Transpose[Conjugate[b]].s.b = I
        Complex bnorm = 0.;

        for(row=0; row<6; row++) {
            for(Subscript col=0; col<6; col++) {
                bnorm += Complex(real(b(row)),-imag(b(row))) * Complex(s(row,col),0) * b(col);
            }
        }
        for(row=0; row<6; row++) {
            b(row) *= sqrt(Complex(1,0)/abs(bnorm));
        }

        if(imag(bnorm)<0) {
            lambda = Complex(real(lambda),-imag(lambda));
            for(int row=0; row<6; row++)
                b(row) = Complex(real(b(row)),-imag(b(row)));
        }

        eigenvalues(pln) = lambda;
        eigenvectors.row(pln) = b;
    }
}

void Symplectify(RealMatrix& a)
{
    IdentityMatrix Ident(6);
    RealMatrix I(Ident);
    RealMatrix s(Ident);

    for(int n=0; n<3; n++)
    {
        int m = 2 * n;
        s(m,m) = s(m+1,m+1) = 0;
        s(m,m+1) = 1;
        s(m+1,m) = -1;
    }

    RealMatrix Ipa(I + a);
#ifndef NDEBUG
    MatrixForm(Ipa,debug_os);
#endif
    Inverse(Ipa);

    RealMatrix Ima(I - a);
    RealMatrix v(s * Ima * Ipa);

    RealMatrix w(6);
    for(int row=0; row<6; row++)
        for(int col=0; col<6; col++)
            w(row,col) = (v(row,col) + v(col,row))/2;

    RealMatrix sw(s * w);

    a = I - sw;
    Inverse(a);
    a = a * (I + sw);

    //Check Symplecticity
    /*	RealMatrix aT(6);
    for(row=0; row<6; row++)
    for(int col=0; col<6; col++)
    aT(row,col) = a(col,row);

      RealMatrix aTsa(aT * s * a);
      ofstream smatrix("smatrix.dat");
      for(row=0; row<6; row++)
      {
      for(int col=0; col<6; col++) smatrix<<std::setw(14)<<aTsa(row,col);
      smatrix<<endl;
      }
    */
}

void tred2(RealMatrix& a, RealVector& d, RealVector& e)
{
    int n = a.nrows();
    int i,j,k;
    d.redim(n);
    e.redim(n);

    for(i=n-1; i>0; i--)
    {
        double h = 0;
        double scale = 0;
        if(i>1)
        {
            for(k=0; k<i; k++) scale += fabs(a(i,k));
            if(scale==0)
                e(i) = a(i,i-1);
            else
            {
                for(k=0; k<i; k++)
                {
                    a(i,k) /= scale;
                    h += a(i,k) * a(i,k);
                }

                double f = a(i,i-1);
                double g = (f>=0.0) ? -sqrt(h) : sqrt(h);
                e(i) = scale * g;
                h -= f * g;
                a(i,i-1) = f - g;
                f = 0.0;

                for(j=0; j<i; j++)
                {
                    a(j,i) = a(i,j) / h;
                    g = 0.0;

                    for(k=0;  k<=j; k++) g += a(j,k) * a(i,k);
                    for(k=j+1; k<i; k++) g += a(k,j) * a(i,k);

                    e(j) = g / h;

                    f += e(j) * a(i,j);
                }

                double hh = f / (h+h);

                for(j=0; j<i; j++)
                {
                    f = a(i,j);
                    e(j) = g = e(j) - hh * f;

                    for(k=0; k<=j; k++) a(j,k) -= (f * e(k) + g * a(i,k));
                }
            }
        }
        else
            e(i) = a(i,i-1);

        d(i) = h;
    }

    d(0) = 0.0;
    e(0) = 0.0;

    for(i=0; i<n; i++)
    {
        if(d(i))
            for(j=0; j<i; j++)
            {
                double g = 0.0;
                for(int k=0; k<i; k++) g += a(i,k) * a(k,j);
                for(int k=0; k<i; k++) a(k,j) -= g * a(k,i);
            }

        d(i) = a(i,i);
        a(i,i) = 1.0;

        for(j=0; j<i; j++) a(j,i) = a(i,j) = 0.0;
    }
}

double pythag(double a, double b)
{
    double absa = fabs(a);
    double absb = fabs(b);
    if( absa > absb)
        return absa * sqrt(1.0 + pow(absb/absa,2));
    else
        return (absb == 0.0) ? 0.0 : absb * sqrt(1.0 + pow(absa/absb,2));
}

void tqli(RealVector& d, RealVector& e, RealMatrix& z)
{
    int m = 0;
    int n = d.size();
    int i,k,l;

    for(i=1; i<n; i++) e(i-1) = e(i);
    e(n-1) = 0.0;

    for(l=0; l<n; l++)
    {
        int iter = 0;
        do
        {
            for(m=l; m<n-1; m++)
            {
                double dd = fabs(d(m)) + fabs(d(m+1));
                if(!e(m)) break;
            }

            if(m != l)
            {
                double g = (d(l+1) - d(l))/(2.0*e(l));
                double r = pythag(g,1.0);
                g = d(m) - d(l) + e(l) / (g + _SIGN(r,g));
                double s = 1.0;
                double c = 1.0;
                double p = 0.0;

                for(i=m-1; i>=l; i--)
                {
                    double f = s * e(i);
                    double b = c * e(i);
                    r = pythag(f,g);
                    e(i+1) = r;
                    if(r==0.0)
                    {
                        d(i+1) -= p;
                        e(m) = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d(i+1) - p;
                    r = (d(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d(i+1) = g + p;
                    g = c * r - b;

                    for(k=0; k<n; k++)
                    {
                        f = z(k,i+1);
                        z(k,i+1) = s * z(k,i) + c * f;
                        z(k,i) = c * z(k,i) - s * f;
                    }
                }

                if(r==0.0 && i>=l) continue;

                d(l) -= p;
                e(l) = g;
                e(m) = 0.0;
            }
        } while( (m!=l) && (iter++<30) );
    }
}

void EigenSystemSymmetricMatrix(RealMatrix& m, RealVector& eigenvalues)
{
    int n = m.nrows();
    eigenvalues.redim(n);

    RealVector e(n);
    tred2(m, eigenvalues, e);
    tqli(eigenvalues, e, m);
}

}; // end namespace TLAS

