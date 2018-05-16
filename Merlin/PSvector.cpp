/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/22 09:28:00 $
// $Revision: 1.3 $
//
/////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include "PSvector.h"

using namespace std;

ostream& operator<<(ostream& os, const PSvector& v)
{
	for(size_t i=0; i<PS_LENGTH; i++)
	{
		os<<setw(24)<<scientific<<setprecision(10)<<v[i];
	}

//    copy(v.v,v.v+6,ostream_iterator<double>(os," "));
	return os<<'\n';
}

std::istream& operator>>(std::istream& is, PSvector& v)
{
	for(double *q = v.v; (q!=v.v+PS_LENGTH) && is; q++)
	{
		is>>*q;
	}
	return is;
}

