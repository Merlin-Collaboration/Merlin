/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef stdext_algorithm_h
#define stdext_algorithm_h 1

#include "merlin_config.h"

template<class InIt, class OutIt, class Pred>
OutIt copy_if(InIt first, InIt last, OutIt x, Pred test)
{
	while(first++ != last) 
		if(test(*first)) *x++ = *first;
	return x;
}

#endif
