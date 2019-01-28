/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "MatrixPrinter.h"
#include <iomanip>
#include <vector>

using namespace std;

typedef vector<string> StringArray;

namespace
{
inline int max(int a, int b)
{
	return (a > b) ? a : b;
}

} // end anonymous namespace

void MatrixForm(const RealMatrix& M, std::ostream& os, const OPFormat& fmt)
{
	StringArray sarray;
	int width = 0;
	Subscript i, j;

	sarray.reserve(M.nrows() * M.ncols());

	for(i = 0; i < M.nrows(); i++)
	{
		for(j = 0; j < M.ncols(); j++)
		{
			sarray.push_back(fmt(M(i, j)));
			width = max(width, static_cast<int>(sarray.back().size()));
		}
	}

	// single space padding
	width++;

	StringArray::iterator s = sarray.begin();

	for(i = 0; i < M.nrows(); i++)
	{
		for(j = 0; j < M.ncols(); j++, s++)
		{
			os << setw(width) << right << (*s).c_str();
		}
		os << endl;
	}
}
