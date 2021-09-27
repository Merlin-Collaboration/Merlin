/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_TLASimp
#define _h_TLASimp

// Templated functions and classes for performing linear algebra. The following
// perform linear algebra functions based on TMAT::Matrix and TMAT::Vector classes.

#include "TLAS.h"

namespace TLAS
{
// Class Declarations

/**
 * \class LUMatrix
 *
 * Represents a LU-decomposition of a square matrix. An LUMatrix object
 * can be used to "solve" a system of linear equations with many RHS vectors.
 */
template<class T>
class LUMatrix
{
public:

	/**
	 * Construction from an arbitrary (square) matrix.
	 * Throws DimensionError() if the matrix is not square, or
	 * SingularMatrix() if the matrix is singular.
	 */
	explicit LUMatrix(const Matrix<T>& M) :
		lud(M), indexes(), d(1)
	{
		LUfactor(lud, indexes, d);
	}

	/**
	 * Solve the RHS vector using this LU-decomposition.
	 * Note that rhs is over-written with the result.
	 * Returns rhs.
	 */
	Vector<T>& operator()(Vector<T>& rhs) const
	{
		return LUsolve(lud, indexes, rhs);
	}
	SubVector<T>& operator()(SubVector<T>& rhs) const
	{
		return LUsolve(lud, indexes, rhs);
	}

	/**
	 * These functions are provided for acting on const. rhs vectors.
	 */
	Vector<T> operator()(const Vector<T>& rhs) const
	{
		Vector<T> v(rhs);
		return operator()(v);
	}
	Vector<T> operator()(const SubVector<T>& rhs) const
	{
		Vector<T> v(rhs);
		return operator()(v);
	}
	Vector<T> operator()(const ConstSubVector<T>& rhs) const
	{
		Vector<T> v(rhs);
		return operator()(v);
	}

	/**
	 * Return the determinant of the original matrix
	 */
	T det() const
	{
		double dt = d;
		for(int i = 0; i < lud.ncols(); i++)
		{
			dt *= lud(i, i);
		}
		return dt;
	}

private:
	Matrix<T> lud;
	std::vector<int> indexes;
	T d;
};

/**
 * \class SVDMatrix
 *
 * Represents a Singular-Value Decomposition of a matrix. As with LUMatrix, an
 * SVDMatrix object can be used to solve a system of linear equations with many
 * RHS vectors.
 */

template<class T>
class SVDMatrix
{
public:

	/**
	 * Construction from an arbitrary matrix. If M.nrows()<M.ncols(), then M is first
	 * made square by the addition of M.ncols()-M.nrows() zero rows. threshold specifies
	 * the relative (to the largest singular value) threshold value below which the
	 * singular values are set to zero.
	 */
	explicit SVDMatrix(const Matrix<T>& M, T threshold = T(1e-06));
	SVDMatrix(const Matrix<T>& M, const Vector<T>& wts, T threshold = T(1e-06));

	/**
	 * Solve the RHS vector using this SVD.
	 */
	Vector<T> operator()(const Vector<T>& rhs) const
	{
		Vector<T> x(w.size());
		Vector<T> y(rhs);
		y *= wts;
		return SVDsolve(u, w, v, y, x);
	}
	Vector<T> operator()(const SubVector<T>& rhs) const
	{
		Vector<T> x(w.size());
		Vector<T> y(rhs);
		y *= wts;
		return SVDsolve(u, w, v, y, x);
	}
	Vector<T> operator()(const ConstSubVector<T>& rhs) const
	{
		Vector<T> x(w.size());
		Vector<T> y(rhs);
		y *= wts;
		return SVDsolve(u, w, v, y, x);
	}
	const std::vector<bool>& wflags() const
	{
		return wflgs;
	}

	/**
	 * Decomposed matrices
	 */
	const Matrix<T>& U() const
	{
		return u;
	}
	const Matrix<T>& V() const
	{
		return v;
	}
	const Vector<T>& W() const
	{
		return w;
	}

private:

	void Init(const Matrix<T>& M, T threshold);

	Matrix<T> u, v;
	Vector<T> w;
	Vector<T> wts;
	std::vector<bool> wflgs;
};

template<class T>
SVDMatrix<T>::SVDMatrix(const Matrix<T>& M, T threshold) :
	wts(M.nrows())
{
	wts = T(1);
	Init(M, threshold);
}

template<class T>
SVDMatrix<T>::SVDMatrix(const Matrix<T>& M, const Vector<T>& wts1, T threshold) :
	wts(M.nrows())
{
	wts = wts1;
	Init(M, threshold);
}

template<class T>
void SVDMatrix<T>::Init(const Matrix<T>& M, T threshold)
{
	if(M.nrows() < M.ncols())
	{
		u.redim(M.ncols(), M.ncols());
		u(Range(0, M.nrows() - 1), Range(0, M.ncols() - 1)) = M;
	}
	else
	{
		u.copy(M);
	}
        

	// adjust for weights: multiply each column vector by wts
	for(Subscript ir = 0; ir < u.ncols(); ir++)
	{
		u.column(ir) *= wts;
	}

	w.redim(u.ncols());
	wflgs = std::vector<bool>(u.ncols(), true);

	v.redim(u.ncols(), u.ncols());
        SVD(u, w, v);

	T wmin = threshold != T(0) ? threshold * (*std::max_element(w.begin(), w.end())) : threshold;
	int zerocount = 0;
	for(size_t i = 0; i < w.size(); i++)
		if(w[i] <= wmin)
		{
			w[i] = T(0);
			wflgs[i] = false;
			zerocount++;
		}
	if(zerocount == static_cast<int>(w.size()))
	{
		throw SingularValuesAllZero();
	}
}

// Low level routines 

template<class T>
void LUfactor(Matrix<T>& M, std::vector<int>& pivot, T& sign)
{
	int imax = -1;

	if(M.ncols()  != M.nrows())
	{
		throw DimensionError();
	}
	const int n = M.ncols();

	std::vector<T> scale(n);

	sign = 1.0;
	for(int i = 0; i < n; i++)
	{
		double big = 0.0;
		for(int j = 0; j < n; j++)
		{
			big = std::max(fabs(M[i][j]), big);
		}
		if(big == 0.0)
		{
			throw SingularMatrix();
		}
		scale[i] = 1.0 / big;
	}
	for(int j = 0; j < n; j++)
	{
		for(int i = 0; i < j; i++)
		{
			for(int k = 0; k < i; k++)
			{
                                M[i][j]-= M[i][k]*M[k][j]; 
			}
		}
		double big = 0.0;
		for(int i = j; i < n; i++)
		{
			for(int k = 0; k < j; k++)
			{
                                M[i][j]-= M[i][k]*M[k][j]; 
			}
                        double howbig = scale[i] * fabs(M[i][j]);
			if( howbig >= big)
			{
				big = howbig;
				imax = i;
			}
		}
		if(j != imax)
		{
			for(int k = 0; k < n; k++)
			{
				double dum = M[imax][k];
				M[imax][k] = M[j][k];
				M[j][k] = dum;
			}
			sign = - sign; 
			scale[imax] = scale[j];
		}
		pivot[j] = imax;
		if(M[j][j] == 0.0)
		{
			M[j][j] = std::numeric_limits<T>::epsilon();
		}
		if(j < n)
		{
			double dum = 1.0 / (M[j][j]);
			for(int i = j + 1; i < n; i++)
			{
				M[i][j] *= dum;
			}
		}
	}
}

// LU back substitution

template<class T, class V>
V& LUsolve(const Matrix<T>& M, const std::vector<int>& pivot, V& v)
{
/*   
 *   M  and pivot give the LU decomposed matrix
 *   v is the RHS in the equation Mx=v   overwritten by the solution
 *   */
	int  first = -1; // first non-vanishing element of v
	const int n = M.ncols();
	T sum;
// Do forward substitution
//
	for(int i = 0; i < n; i++)
	{
		sum = v[pivot[i]];
		v[pivot[i]] = v[i];
		if(first != -1)
		{
			for(int j = first; j < i ; j++)
			{
				sum -= M[i][j] * v[j];
			}
		}
		else if(sum != 0)
		{
			first = i;
		}
		v[i] = sum;
	}
   
// Now do back substitution
//
	for(int i = n - 1; i >= 0; i--)
	{
		sum = v[i];
		for(int j = i + 1; j < n; j++)
		{
			sum -= M[i][j] * v[j];
		}
		v[i] = sum / M[i][i];
	}
	return v;
}

template<class T>
Matrix<T>& InvertMatrix(Matrix<T>& m)
{
	if(m.nrows() != m.ncols())
	{
		throw NonSquareMatrix();
	}

	LUMatrix<T> lu(m);
	for(int i = 0; i < m.ncols(); i++)
	{
		SubVector<T>& col = m.column(i);
		col = T(0);
		col[i] = T(1);
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



template<class T, class V>
Vector<T>& SVDsolve(const Matrix<T>& u, const Vector<T>& w, const Matrix<T>& v, const V& b, Vector<T>& x)
{
     Vector<T> work=Transpose(u) * b;
     for(uint i=0;i<work.size();i++) work[i]= (w[i]==0.0) ? 0 : work[i]/w[i]; // need to catch divide by zero
     x = v * work;
    return x;
}

template<class T>
inline T PYTHAG(T a, T b)
{
	a = fabs(a);
	b = fabs(b);
	if(a > b)
	{
		T c = b / a;
		return a * sqrt(1 + c * c);
	}
	else if(!fequal(b, 0.0))
	{
		T c = a / b;
		return b * sqrt(1 + c * c);
	}
	else
	{
		return T(0);
	}
}

template<class T> inline T SIGN(T a, T b)
{
	return b >= T(0) ? fabs(a) : -fabs(a);
}

template<class T>
void SVD(Matrix<T>& M, Vector<T>& V, Matrix<T>& S)
{
	assert(M.nrows() <= M.ncols());
        const int rmax = M.nrows() -1 ;
        const int cmax = M.ncols() -1 ;

	T f0, f1, f2, f3;
	T F0 = 0.0, F1 = 0.0, F2 = 0.0;

	Vector<T> work(cmax+1);

	for(int i = 0; i <= cmax ; i++)
	{
                Range R1(i,rmax);
		work[i] = F0 * F1;
		F1 = f1 = F0 = 0.0;
		if(i <= rmax)
		{
			F0 += M(R1,i).sumfabs();
			if(!fequal(F0, 0.0))
			{
                                M(R1,i) /= F0;
                                f1 += M(R1,i).sumsquared() ;
				f0 = M[i][i];
				F1 = -SIGN(sqrt(f1), f0);
				f3 = f0 * F1 - f1;
				M[i][i] = f0 - F1;
				if(i != cmax)
				{
					for(int j = i+1; j <= cmax; j++)
					{
                            
						f1 = M(R1,i) * M(R1,j);
						f0 = f1  / f3;
                                                M(R1,j) +=  f0*M(R1,i);
					}
				}
                                M(R1,i) *= F0;
			}
		}
		V[i] = F0 * F1;
		F1 = f1 = F0 = 0.0;
		if(i <= rmax  && i != cmax )
		{
                        Range R2(i+1,cmax);
                        F0 += M(i,R2).sumfabs();
			if(!fequal(F0, 0.0))
			{
                                M(i,R2) /= F0;
                                f1 += M(i,R2).sumsquared() ;
				f0 = M[i][i+1];
				F1 = -SIGN(sqrt(f1), f0);
				f3 = f0 * F1 - f1;
				M[i][i+1] = f0 - F1;
				work(R2) = 0; // assignment to subvector from subvector is private for some reason 
                                work(R2) += M(i,R2) ;
                                work(R2) /=f3;
				if(i != rmax )
				{
					for(int j = i+1; j <= rmax; j++)
					{
                                                f1= M(j,R2) * M(i,R2);
				                M(j,R2) +=  f1 * work(R2) ;
					}
				}
				M(i,R2) *= F0;
			}
		}
		F2 = std::max(F2, (fabs(V[i]) + fabs(work[i])));
	}
	for(int i = cmax ; i >= 0; i--)
	{
		if(i < cmax )
		{
                        Range R2(i+1,cmax);
			if(!fequal(F1, 0.0))
			{
				S(R2,i) = (M(i,R2) / M[i][i+1]) /F1;
				for(int j = i+1; j <= cmax ; j++)
				{
					f1 = M(i,R2) * S(R2,j); 
                                        S(R2,j) += f1 * S(R2,i);
				}
			}
                        S(i,R2)=0;
                        S(R2,i)=0;
		}
		S[i][i] = 1.0;
		F1 = work[i];
	}
	for(int i = cmax ; i >= 0; i--)
	{
		Range R3(i+1,rmax); 
		Range R4(i+1,cmax); 
		Range R5(i,rmax); 
		F1 = V[i];
		if(i < cmax )   M(i,R4)=0;
		if(!fequal(F1, 0.0))
		{
			F1 = 1.0 / F1;
			if(i != cmax )
			{
				for(int j = i+1; j <= cmax ; j++)
				{
                                        f1=M(R3,i) * M(R3,j);
					f0 = (f1 / M[i][i]) * F1;
                                        M(R5,j) += f0*M(R5,i);
				}
			}
                        M(R5,i) *= F1;
		}
		else
		{
                        M(R5,i) = 0; 
		}
		++M[i][i];
	}
	for(int k = cmax ; k >= 0; k--)
	{
		for(int ntry = 30; ntry >= 0; ntry--)
		{
			if(ntry == 0)
			{
				throw ConvergenceFailure();
			}
			bool ok=true;
                        int L,remember;
			for( L = k; L >= 0; L--)
			{
				remember = L - 1;
				if(fequal(fabs(work[L]) + F2, F2))
				{
					ok=false; 
					break;
				}
				if(fequal(fabs(V[L-1]) + F2, F2))
				{
					break;
				}
			}
			if(ok) 
			{
				f1 = 1.0;
				f2 = 0.0;
				for(int i = L; i <= k; i++)
				{
					f0 = f1 * work[i];
					if(!fequal(fabs(f0) + F2, F2))
					{
						F1 = V[i];
						f3 = PYTHAG(f0, F1);
						V[i] = f3;
						f3 = 1.0 / f3;
						f2 = F1 * f3;
						f1 = (-f0 * f3);
                                                Vector<T> Y=M.column(remember);
                                                Vector<T> Z=M.column(i);
                                                M.column(remember)= Y * f2 + Z * f1;
                                                M.column(i) = Z * f2 - Y * f1;
					}
				}
			}
			T q1 = V[k];
			if(L == k)
			{
				if(q1 < 0.0)
				{
					V[k] = -q1;
                                        S(Range(0,cmax),k) *= -1;
				}
				break;
			}

			T q2 = V[L];
			T q3 = V[k-1];
			F1 = work[k-1];
			f3 = work[k];
			f0 = ((q3 - q1) * (q3 + q1) + (F1 - f3) * (F1 + f3)) / (2.0 * f3 * q3);
			F1 = PYTHAG(f0, 1.0);
			f0 = ((q2 - q1) * (q2 + q1) + f3 * ((q3 / (f0 + SIGN(F1, f0))) - f3)) / q2;
			f2 = f1 = 1.0;
			for(int j = L; j <= (k-1); j++)
			{
				F1 = work[j+1];
				q3 = V[j+1];
				f3 = f1 * F1;
				F1 = f2 * F1;
				q1 = PYTHAG(f0, f3);
				work[j] = q1;
				f2 = f0 / q1;
				f1 = f3 / q1;
				f0 = q2 * f2 + F1 * f1;
				F1 = F1 * f2 - q2 * f1;
				f3 = q3 * f1;
				q3 = q3 * f2;
                                Vector<T> X=S.column(j);
                                Vector<T> Z1=S.column(j+1);
                                S.column(j) = X * f2 + f1 * Z1;
                                S.column(j+1) = Z1 * f2 - X * f1;
				q1 = PYTHAG(f0, f3);
				V[j] = q1;
				if(!fequal(q1, 0.0))
				{
					q1 = 1.0 / q1;
					f2 = f0 * q1;
					f1 = f3 * q1;
				}
				f0 = (f2 * F1) + (f1 * q3);
				q2 = (f2 * q3) - (f1 * F1);
                                Vector<T> Y=M.column(j);
                                Vector<T> Z=M.column(j+1);
                                M.column(j)= f2 * Y + f1 * Z;
                                M.column(j+1)= f2 * Z - f1 * Y;
			}
			work[L] = 0.0;
			work[k] = f0;
			V[k] = q2;
		}
	}
}

} // end namespace TLAS

#endif // _h_TLASimp
