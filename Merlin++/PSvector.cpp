/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iomanip>
#include "PSvector.h"

using namespace std;

ostream& operator<<(ostream& os, const PSvector& v)
{
	for(size_t i = 0; i < PS_LENGTH; i++)
	{
		os << setw(24) << scientific << setprecision(10) << v[i];
	}

//    copy(v.v,v.v+6,ostream_iterator<double>(os," "));
	return os << '\n';
}

std::istream& operator>>(std::istream& is, PSvector& v)
{
	for(double *q = v.v; (q != v.v + PS_LENGTH) && is; q++)
	{
		is >> *q;
	}
	return is;
}
