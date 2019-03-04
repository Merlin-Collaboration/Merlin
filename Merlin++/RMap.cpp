/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "RMap.h"
#include "MatrixPrinter.h"
#include "LinearAlgebra.h"
#include <algorithm>

void RMap::ToMatrix(RealMatrix& R, bool init) const
{
	if(init)
	{
		R = RealMatrix(6, 6, 0);
	}

	for(const_itor t = rterms.begin(); t != rterms.end(); t++)
	{
		R(t->i, t->j) = t->val;
	}
}

RMap::operator RealMatrix() const
{
	RealMatrix R(6, 6, 0);
	ToMatrix(R, false);
	return R;
}

void RMap::MatrixForm(std::ostream& os) const
{
	RealMatrix R(6, 6, 0.0);
	ToMatrix(R, false);
	::MatrixForm(R, os);
}
/*
   PSvector& RMap::Apply(PSvector& X) const
   {
    PSvector Y(0);
    Apply(X,Y);
    return X=Y;
   }
 */
void RMap::Apply(const PSvector& orig, PSvector& result) const
{
	for(const_itor r = rterms.begin(); r != rterms.end(); r++)
	{
		r->Apply(orig, result);
	}
}

RMap::itor RMap::FindTerm(int i, int j)
{
	itor ri = rterms.begin();
	for(; ri != rterms.end(); ri++)
		if(ri->i == i && ri->j == j)
		{
			break;
		}
	return ri;
}

double RMap::operator()(int i, int j) const
{
	const_itor ri = const_cast<RMap*>(this)->FindTerm(i - 1, j - 1);
	return ri != rterms.end() ? ri->val : 0;
}

double& RMap::operator()(int i, int j)
{
	i--;
	j--;
	itor ri = FindTerm(i, j);
	if(ri == rterms.end())
	{
		rterms.push_back(Rij(i, j, 0));
		ri = rterms.end();
		ri--;
	}
	return ri->val;
}

RMap::RMap(const RealMatrix& R) :
	rterms()
{
	rterms.reserve(8);
	for(size_t i = 0; i < R.nrows(); i++)
		for(size_t j = 0; j < R.ncols(); j++)
		{
			if(R(i, j) != 0)
			{
				rterms.push_back(Rij(i, j, R(i, j)));
			}
		}
}

void MakeIdentity(RMap& R, int ndf)
{
	for(int i = 1; i <= 2 * ndf; i++)
	{
		R.AddTerm(i, i, 1.0);
	}
}
