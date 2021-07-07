/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _MultiNormal_h
#define _MultiNormal_h

#include "MerlinException.h"
#include "LinearAlgebra.h"
#include "TCovMtrx.h"
#include "RandomNG.h"

#include "MatrixPrinter.h"
#include "OPFormat.h"

// DK 2006-02-23

/**
 * A multivariate normal random number generator
 * by Cholesky decomposition: L*Lt=Cov, v=L*x, x<-normal(0,1)
 * input: covariance and mean
 * output: N-dim normal random number according to covariance
 * RandomNG must be initialized
 */

template<int N>
class MultiNormal
{
public:
	MultiNormal(const RealMatrix& Cov, const RealVector& M);
	MultiNormal(const TCovMtrx<double, N>& Cov, const RealVector& M);
	MultiNormal(const TPSMoments<N / 2>& pm);
	~MultiNormal();
	RealVector GetRandVec();   // get a random vector v=Lx+m
private:
	void  CholeskyDecomp();    // do Cov->L
	const RealVector Mean;    // vector(N) of mean values
	RealMatrix L;    // matrix(NxN): covariance in upper tri. filled by constructor
	// cholesky decomp of covariance matrix in lower tri. filled by CholeskyMatrix
	double*              x;    // pointer to scratch pad for random numbers

	//Copy protection
	MultiNormal(const MultiNormal& rhs);
	MultiNormal& operator=(const MultiNormal& rhs);
};

template<int N>
inline RealVector MultiNormal<N>::GetRandVec()    // v=Lx+m
{
	for(int i = 0; i < N; i++)
	{
		x[i] = RandomNG::normal(0, 1);
	}
	for(int i = N - 1; i >= 0; i--)
	{
		x[i] *= L(i, i);
		x[i] += Mean(i);
		for(int j = 0; j < i; j++)
		{
			x[i] += L(i, j) * x[j];
		}
	}
	return RealVector(x, N);
}

template<int N>
MultiNormal<N>::MultiNormal(const RealMatrix& Cov, const RealVector& M) :
	Mean(M), L(Cov)
{
	if(N == 0 || N > Cov.ncols() || Cov.ncols() != Cov.nrows() || N > Mean.size())
	{
		throw DimensionError();
	}
	CholeskyDecomp();
	// MatrixForm(L,cout,OPFormat().precision(6).fixed());
	x = new double[N];
}

template<int N>
MultiNormal<N>::MultiNormal(const TCovMtrx<double, N>& Cov, const RealVector& M) :
	Mean(M), L(N, N)
{
	if(N == 0 || N > Mean.size())
	{
		throw DimensionError();
	}
	for(int i = 0; i < N; i++)
		for(int j = i; j < N; j++)
		{
			L(i, j) = Cov(i, j);
		}
	CholeskyDecomp();
	// MatrixForm(L,cout,OPFormat().precision(6).fixed());
	x = new double[N];
}

template<int N>
MultiNormal<N>::MultiNormal(const TPSMoments<N / 2>& pm) :
	Mean(pm), L(N, N)
{
	if(N == 0 || N > Mean.size())
	{
		throw DimensionError();
	}
	for(int i = 0; i < N; i++)
	{
		for(int j = i; j < N; j++)
		{
			L(i, j) = pm.sig(i, j);
		}
	}
	CholeskyDecomp();
	x = new double[N];
}

template<int N>
MultiNormal<N>::~MultiNormal()
{
	delete[] x;
}

template<int N>
void MultiNormal<N>::CholeskyDecomp()// uses Cholesky-Banachiewicz algorithm from Wikipedia article 
	//   en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm
	//   The original symmetric matrix is given in the upper right part of L
	//   The desired matrix is put in the lower left part (including the diagonal)
{
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			double stot=0;
			for(int k=0;k<j;k++)
			{
				stot += L(i,k)*L(j,k);
			}
			if (i==j) 
			{
				double x=L(i,j)-stot;
				if(x<0) throw MerlinException("Negative square root in CholeskyDecomp");
				L(i,j)=sqrt(x);

			} 
			else
			{
				L(i,j)=(L(j,i)-stot)/L(j,j);
			}
		}
	}
}

#endif
