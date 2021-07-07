/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "utils.h"
#include "NumericalConstants.h"

double LogGamma(double xx)
{
	//   Use Lanczos approximation

	// code taken from Wikipedia    en.wikipedia.oirg/wiki/Lanczos_approximation
	// 
	const int N=8;
	static const double p[N]={676.5203681218851,-1259.1392167224028, 771.32342877765313,-176.61502916214059, 12.507343278686905,
		-0.13857109526572012, 9.9843695780195716E-06, 1.5056327351493116E-7};
	const double c=log(sqrt(2*pi));

	double z=xx-1;
	double sum=0.99999999999980993;
	for(int i=0; i<N;i++) sum+= p[i]/(z+i+1);
	double val=c+(z+0.5)*log(z+N-.5)-(z+N-.5) + log(sum) ;
	 //std::cout<<" new style " <<val<<std::endl;
	return  val;
}
