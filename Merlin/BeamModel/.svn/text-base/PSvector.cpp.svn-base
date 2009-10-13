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

//#include <iterator>
//#include <algorithm>
#include <iomanip>
#include "BeamModel/PSvector.h"

using namespace std;

ostream& operator<<(ostream& os, const PSvector& v)
{
    for(size_t i=0; i<6; i++)
        os<<setw(14)<<scientific<<setprecision(5)<<v[i];

//    copy(v.v,v.v+6,ostream_iterator<double>(os," "));
    return os<<'\n';
}

std::istream& operator>>(std::istream& is, PSvector& v)
{
    for(double *q = v.v; (q!=v.v+6) && is; q++)
        is>>*q;
    return is;
}

