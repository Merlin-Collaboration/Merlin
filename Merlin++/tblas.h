/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_tblas
#define _h_tblas

#include "utils.h"
#include <cassert>

namespace tblas1
{

// dot <- alpha + x'.y
template<class V1, class V2, class T>
T tdot(const V1& v1, const V2& v2, T alpha)
{
	for(size_t i = 0; i < v1.size(); ++i)
	{
		alpha += v1(i) * v2(i);
	}
	return alpha;
}

template<class V1, class V2, class T>
void taxpy(T alpha, const V1& x, V2& y)
{
	assert(&x != &y);
	if(alpha == 0)
	{
		return;
	}
	if(alpha == 1)
	{
		y += x;
	}
	else
	{
		for(size_t i = 0; i < x.size(); ++i)
		{
			y(i) += alpha * x(i);
		}
	}
}
}

namespace tblas2
{

using namespace tblas1;

// y <- alpha*A.x + beta*y
template<class Ta, class M, class Vx, class Tb, class Vy>
void tgemv(bool t, const Ta& alpha, const M& A, const Vx& x, const Tb& beta, Vy& y)
{
// This seems to be where all the work gets done - JM

	assert(&x != &y); // that's a no-no!
	if(!fequal(beta, 1.0))
	{
		for(size_t i = 0; i < A.nrows(); i++)
		{
			y[i] *= beta;
		}
	}

	if(!fequal(alpha, 0.0))
	{
		//start if/else - should be correct
		if(t)
		{
			// use transpose of A
			for(size_t i = 0; i < A.ncols(); ++i)
			{
				for(size_t j = 0; j < A.nrows(); ++j)
				{
					y[i] += alpha * A(j, i) * x[j];
				}
			}
		}
		else
		{
			// use normal A
			for(size_t i = 0; i < A.nrows(); ++i)
			{
				for(size_t j = 0; j < A.ncols(); ++j)
				{
					y[i] += alpha * A(i, j) * x[j];
				}
			}
		}
	} //end
}

// add Symm, Upper and Lower diagonal forms later

// A <- alpha*x.y' + A (outer product)
template<class T, class Vx, class Vy, class M>
void tger(const T& alpha, const Vx& x, const Vy& y, M& A)
{
	if(alpha == 0)
	{
		return;
	}

	for(size_t i = 0; i < A.nrows(); ++i)
		for(size_t j = 0; j < A.ncols(); ++j)
		{
			A(i, j) += alpha * x(i) * y(j);
		}
}
}

namespace tblas3
{

using namespace tblas2;

// C <- alpha*A.B + beta*C
template<class Ta, class Ma, class Mb, class Tb, class Mc>
void tgemm(bool tpa, bool tpb, const Ta& alpha, const Ma& A, const Mb& B, const Tb& beta, Mc& C)
{
	assert(&C != &A && &C != &B); // that's a no-no

	if(fequal(beta, 0.0))
	{
		C = 0;    // assignment faster than multiplication?
	}
	else if(!fequal(beta, 1.0))
	{
		C *= beta;
	}

	if(fequal(alpha, 0.0))
	{
		return;
	}

	// calculate the range of the k subscript
	const size_t klim = tpa ? A.nrows() : A.ncols();

	for(size_t i = 0; i < C.nrows(); ++i)
		for(size_t j = 0; j < C.ncols(); ++j)
			for(size_t k = 0; k < klim; ++k)
			{
				// 4 cases to deal with here
				if(tpa)
				{
					if(tpb)
					{
						C(i, j) += alpha * A(k, i) * B(j, k);    // A'.B'
					}
					else
					{
						C(i, j) += alpha * A(k, i) * B(k, j);    // A'.B
					}
				}
				else
				{
					if(tpb)
					{
						C(i, j) += alpha * A(i, k) * B(j, k);    // A.B'
					}
					else
					{
						C(i, j) += alpha * A(i, k) * B(k, j);    // A.B
					}
				}
			}
}

// Matrix rotation C <- R.M.R'
template<class Tr, class Ta, class Tb>
void tgemr(const Tr& R, const Ta& M, Tb& C)
{
	const size_t n = R.nrows();
	for(size_t i = 0; i < n; i++)
		for(size_t j = 0; j < n; j++)
		{
			C(i, j) = 0;
			for(size_t k = 0; k < n; k++)
				for(size_t l = 0; l < n; l++)
				{
					C(i, j) += R(i, k) * R(j, l) * M(k, l);
				}
		}
}
}

#endif
