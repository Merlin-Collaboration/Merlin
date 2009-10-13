/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "NumericalUtils/Histogram.h"

using std::vector;

size_t Hist(const vector<double>& data, double x1, double x2, double dx, vector<double>& hist)
{
    size_t nbins = size_t((x2-x1)/dx)+1;
    hist = vector<double>(nbins,0.0);
    size_t lost=0;

    for(size_t n=0; n<data.size(); n++) {
        if(data[n]<x1)
            lost++;
        else {
            size_t i = (data[n]-x1)/dx;
            if(i<nbins)
                hist[i] += 1.0;
            else
                lost++;
        }
    }
    return lost;
}
