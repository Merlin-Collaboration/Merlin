/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <fstream>
#include "LinearAlgebra.h"
#include "MerlinException.h"
#include <algorithm> // for std::swap
#include <cmath>

#ifndef DEBUG
#include "MatrixPrinter.h"
#endif

using namespace std;

namespace TLAS
{
template class Vector<Complex>;
template class Matrix<Complex>;
template class SubVector<Complex>;
template class SubMatrix<Complex>;
template class ConstSubVector<Complex>;
template class ConstSubMatrix<Complex>;
}

namespace
{

//ofstream debug_os("tlas_debug.dat");

// MSVC++ bug! Compiler should make correct resolution for abs()!!
inline double ABS(double x)
{
	return fabs(x);
}
inline double ABS(const Complex& z)
{
	return abs(z);
}

#define _SIGN(a, b) ((b > 0) ? fabs(a) : -fabs(a))

template<class T> double Inverse(Matrix<T>& a)
{

	if(a.nrows() != a.ncols())
	{
		throw TLAS::NonSquareMatrix();
	}

	int n = a.nrows();
        vector<bool> done(n,false);
        Matrix<T> b = Matrix<T>(IdentityMatrix(n));
	int pcol=0, prow=0;  // Does not need initialising, but avoids compiler warnings
	double minpiv = 1.0e9;

	for(int i = 0; i < n; i++)
	{
		double largest  = 0.0;
		for(int j = 0; j < n; j++)
			if(!done[j]) 
				for(int k = 0; k < n; k++)
				{
				       if(!done[k]) 
					{
						if(ABS(a(j, k)) >= largest)
						{
							largest  = ABS(a(j, k));
							prow = j;
							pcol = k;
						}
					}
				}
                done[pcol]=true;

		if(prow != pcol)
			for(int j = 0; j < n; j++)
			{
				swap(a(prow, j), a(pcol, j));
				swap(b(prow, j), b(pcol, j));
			}

                minpiv = min(minpiv, ABS(a(pcol,pcol)));
		if(minpiv == 0.0)
		{
			throw TLAS::SingularMatrix();
		}

		T pivinv = 1.0 / a(pcol, pcol);
		a(pcol, pcol) = 1.0;
                a.row(pcol) *= pivinv;
                b.row(pcol) *= pivinv;

		for(int j = 0; j < n; j++)
			if(j != pcol)
			{
				T dum = a(j, pcol);
				a(j, pcol) = 0.0;
                                a.row(j) -= dum*a.row(pcol);
                                b.row(j) -= dum*b.row(pcol);
			}
	}
        a = b; 
	return minpiv;
} // end of Inverse function
} // end anonymous namespace

namespace TLAS
{

double Invert(RealMatrix& t)
{
	return Inverse(t);
}

bool EigenSystem(RealMatrix& t, ComplexVector& eigenvalues, ComplexMatrix& eigenvectors)
{
	bool allOK = true;
	ComplexMatrix m(t);
	Matrix<int> s(6, 6, 0);

	s(0, 1) = s(2, 3) = s(4, 5) = 1;
	s(1, 0) = s(3, 2) = s(5, 4) = -1;

	eigenvalues.redim(3);
	eigenvectors.redim(3, 6);

	for(int pln = 0; pln < 3; pln++)
	{
		//Guess an eigenvalue
		Complex cos_mu((t(2 * pln, 2 * pln) + t(2 * pln + 1, 2 * pln + 1)) / 2.0, 0);
		Complex isin_mu = sqrt(cos_mu * cos_mu - Complex(1, 0));
		Complex lambda(cos_mu + isin_mu);

		//Guess an eigenvector
		ComplexVector b(0.0, 6);
		b(2 * pln) = Complex(1 / sqrt(2.0));
		b(2 * pln + 1) = b(2 * pln);

		// Did you accidentally guess right? - added RJB
		ComplexVector Check(0., 6);
		for(int i = 0; i < 6; i++)
			for(int j = 0; j < 6; j++)
				Check(i) += m(i, j) * b(j);
		double test = 0;
		for(int i = 0; i < 6; i++)
			test += pow(ABS(Check(i) - b(i)), 2);
		Subscript n = 0;
		Subscript row = 0;
		if(test > 1.0E-15)           // if the guess is spot on, iteration fails!

		{ //Iterate!
			ComplexMatrix mp(m);

			double prox = 1.0;
			int iter = 0;
			const int maxiter = 50;

			while(prox > 1.0e-16 && iter < maxiter)
			{
				for(n = 0; n < 6; n++)
				{
					mp(n, n) = m(n, n) - lambda;
				}

				ComplexMatrix minv(mp);
				prox = Inverse(minv);

				ComplexVector y(6);
				for(row = 0; row < 6; row++)
				{
					Complex sum = 0.0;
					for(Subscript col = 0; col < 6; col++)
					{
						sum += minv(row, col) * b(col);
					}
					y(row) = sum;
				}

				Complex ynorm = 0.0;
				Complex ydotb = 0.0;
				for(n = 0; n < 6; n++)
				{
					ynorm += y(n) * y(n);
					ydotb += y(n) * b(n);
				}

				b = y / sqrt(ynorm);
				lambda += 1.0 / ydotb;
				iter++;
			}
			if(iter == maxiter)
			{
				cout << "EigenSystem convergence only reached " << prox << endl;
			}
		}          // end of test check

		//Now normalise to Transpose[Conjugate[b]].s.b = I
		Complex bnorm = 0.;

		for(row = 0; row < 6; row++)
		{
			for(Subscript col = 0; col < 6; col++)
			{
				bnorm += Complex(real(b(row)), -imag(b(row))) * Complex(s(row, col), 0) * b(col);
			}
		}

		for(row = 0; row < 6; row++)
		{
			b(row) *= sqrt(Complex(1, 0) / abs(bnorm));
		}

		if(imag(bnorm) < 0)
		{
			lambda = Complex(real(lambda), -imag(lambda));
			for(int row = 0; row < 6; row++)
			{
				b(row) = Complex(real(b(row)), -imag(b(row)));
			}
		}
		eigenvalues(pln) = lambda;
		eigenvectors.row(pln) = b;
		if((imag(lambda) == 0) && pln == 2)
		{
			cout << " Unstable in 3rd dimension\n";
			allOK = false;
		}
	}
	return allOK;
}

void Symplectify(RealMatrix& a)
{
	IdentityMatrix Ident(6);
	RealMatrix I(Ident);
	RealMatrix s(Ident);

	for(int n = 0; n < 3; n++)
	{
		int m = 2 * n;
		s(m, m) = s(m + 1, m + 1) = 0;
		s(m, m + 1) = 1;
		s(m + 1, m) = -1;
	}

	RealMatrix Ipa(I + a);

	Inverse(Ipa);

	RealMatrix Ima(I - a);
	RealMatrix v(s * Ima * Ipa);

	RealMatrix w(6);
	for(int row = 0; row < 6; row++)
		for(int col = 0; col < 6; col++)
		{
			w(row, col) = (v(row, col) + v(col, row)) / 2;
		}

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
//  Adapted from the  EISPACK subroutine given below
//
/*
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end

*/

	int n = a.nrows();
	int i, j, k;
	d.redim(n);
	e.redim(n);

	for(i = n - 1; i > 0; i--)
	{
		double h = 0;
		if(i > 1)
		{
		        double scale=a(i,Range(0,i-1)).sumfabs();
			if(scale == 0)
			{
				e(i) = a(i, i - 1);
			}
			else
			{       
                                Range r=Range(0,i-1);
                                a(i,r) /= scale;
                                h +=  a(i,r) * a(i,r) ; 

				double f = a(i, i - 1);
				double g = -copysign(sqrt(h),f); 
				e(i) = scale * g;
				h -= f * g;
				a(i, i - 1) = f - g;
				f = 0.0;

				for(j = 0; j < i; j++)
				{
					a(j, i) = a(i, j) / h;
                                        g = a(j,Range(0,j)) * a(i,Range(0,j)) + a(Range(j+1,i-1),j) * a(i,Range(j+1,i-1));;
					e(j) = g / h;
					f += e(j) * a(i, j);
				}

				double hh = f / (h + h);

				for(j = 0; j < i; j++)
				{
					f = a(i, j);
					e(j) = g = e(j) - hh * f;

					for(k = 0; k <= j; k++)
					{
						a(j, k) -= (f * e(k) + g * a(i, k));
					}
				}
			}
		}
		else
		{
			e(i) = a(i, i - 1);
		}

		d(i) = h;
	}

	d(0) = 0.0;
	e(0) = 0.0;

	for(i = 0; i < n; i++)
	{
		if(d(i))
			for(j = 0; j < i; j++)
			{
                                if(i>0) {
                                      Range r=Range(0,i-1);
                                      double g=a(i,r) * a(r,j);
                                      a(r,j) -= g*a(r,i);
                                }
			}

		d(i) = a(i, i);
		a(i, i) = 1.0;

                if(i>0) {
                      a(i,Range(0,i-1))=0;
                      a(Range(0,i-1),i)=0;
                }
	}
}

double pythag(double a, double b)
{
	double absa = fabs(a);
	double absb = fabs(b);
	if(absa > absb)
	{
		return absa * sqrt(1.0 + pow(absb / absa, 2));
	}
	else
	{
		return (absb == 0.0) ? 0.0 : absb * sqrt(1.0 + pow(absa / absb, 2));
	}
}

void tqli(RealVector& diag, RealVector& off, RealMatrix& M)
{
/*
 * Find the eigenvalues and eigenvectors of a matrix that has been reduced to tridiagonal form
 * diag is the diagonal (from 0 to n-1), off is the off-diagonal (from 1 to n-1) 
 *
 *    adapted from the original EISPACK code copied below
 *
 *
      SUBROUTINE TQLI(D,E,N,NP,Z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
 */
	int n = diag.size();

 // move off-diagonal elements down 1, for convenience.
 //
        for(int i = 1; i < n; i++)
	{
		off(i - 1) = off(i);
	}
	off(n - 1) = 0.0;

	for(int j = 0; j < n; j++)
	{
		int ntries = 0;
	        int m = 0;
		do
		{
			for(m = j; m < n - 1; m++)
			{
                                // testing needs rewriting for double precision
                                //
                                double dd = fabs(diag(m)) + fabs(diag(m+1));
                                double small = 1.E-12 * dd;
                                if(fabs(off(m)) <= small) break;
			}
			if(m != j)
			{
                                if(ntries++ > 30 ) throw MerlinException("Too many iterations in tqli");
				double g = (diag(j + 1) - diag(j)) / (2.0 * off(j));
				double r = pythag(g, 1.0);
				g = diag(m) - diag(j) + off(j) / (g + _SIGN(r, g));
				double s = 1.0;
				double c = 1.0;
				double p = 0.0;

				int i;   // needed outside loop due to break
                                for(i = m - 1; i >= j; i--)
				{
					double f = s * off(i);
					double b = c * off(i);
					r = pythag(f, g);
					off(i + 1) = r;
					if(r == 0.0)
					{
						diag(i + 1) -= p;
						off(m) = 0.0;
						break;
					}

					s = f / r;
					c = g / r;
					g = diag(i + 1) - p;
					r = (diag(i) - g) * s + 2.0 * c * b;
					p = s * r;
					diag(i + 1) = g + p;
					g = c * r - b;

					for(int k = 0; k < n; k++)
					{
						f = M(k, i + 1);
						M(k, i + 1) = s * M(k, i) + c * f;
						M(k, i) = c * M(k, i) - s * f;
					}
				}

				if(r == 0.0 && i >= j)
				{
					continue;
				}

				diag(j) -= p;
				off(j) = g;
				off(m) = 0.0;
			}
		} while(m != j);
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

} // end namespace TLAS
