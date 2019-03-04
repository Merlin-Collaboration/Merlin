/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _MultiNormal_h
#define _MultiNormal_h

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
void MultiNormal<N>::CholeskyDecomp()   // copied from John Ellithorpe based on Numerical Recipeces
{
	double dp[N];
	long i, j, k;
	double sum;
	for(i = 0; i < N; i++)
	{
		for(j = i; j < N; j++)
		{
			sum = L(i, j);
			k = i;
			while(--k >= 0)
			{
				sum -= L(i, k) * L(j, k);
			}
			if(i == j)
			{
				if(sum <= 0.0)
				{
					cout << "MultiNormal::CholeskyDecomp Failed!" << endl;
					// we set L=0 so that we just return the mean values
					for(int i = 0; i < N; i++)
						for(int j = 0; j < N; j++)
						{
							L(i, j) = 0;
						}
					return;
				}
				else
				{
					dp[i] = sqrt(sum);
				}
			}
			else
			{
				L(j, i) = sum / dp[i];
			}
		}
	}
	for(i = 0; i < N; i++)
	{
		L(i, i) = dp[i];
	}
}

#endif
