/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3.2 (2008)
// 
// Copyright: see Merlin/copyright.txt
//
// Created: DK 25.2.2008
//    see BeamDynamics/ParticleTracking/SpoilerWakeProcess.cpp
/////////////////////////////////////////////////////////////////////////


#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/SpoilerPotentialModels.h"
#include <math.h>

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;

//----------------------------------------------------------------------------------------------------------
//   The geometric wake potential
//----------------------------------------------------------------------------------------------------------
TaperedCollimatorPotentials::TaperedCollimatorPotentials(int m, double aa, double bb)
: a(aa), b(bb), SpoilerWakePotentials(m) {

	coeff=new double[m+1];
   	for(int i=0;i<(m+1);i++) {
		coeff[i]=2*(1./pow(a,2*i)-1./pow(b,2*i)); 
	}      

};
TaperedCollimatorPotentials::~TaperedCollimatorPotentials(){
    	if(coeff!=0) delete[] coeff;
};
double TaperedCollimatorPotentials::Wlong (double z, int m) const {
    	return z>0? -((m)/a*coeff[m]/exp((m)*z/a)):0 ;
};
double TaperedCollimatorPotentials::Wtrans(double z, int m) const {
    	return z>0? coeff[m]/exp((m)*z/a):0;
};
//------------------------------------------------------------------------------------------------------------
//      the resistive wake potentials  (in MKS ssytem)
//------------------------------------------------------------------------------------------------------------
ResistiveWakePotentials::ResistiveWakePotentials(int m, double r, double s, double l)
: rad(r), sigma(s), length(l), SpoilerWakePotentials(m) {

	coeff = new double[m+1];

	int delta = 0;
	if (m==0) delta = 1;
	for (int i=0; i<(m+1); i++) {
		coeff[i]= 1/pi*pow(rad,2*i+1)*(1+delta);
	}

};
ResistiveWakePotentials::~ResistiveWakePotentials(){
    	if(coeff!=0) delete[] coeff;
};
double ResistiveWakePotentials::Wlong (double z, int m) const {
   	return z>0? coeff[m]*sqrt(1/sigma*376.6)*sqrt(z)*length:0;
};
double ResistiveWakePotentials::Wtrans(double z, int m) const {
   	return z>0? -2*coeff[m]*sqrt(SpeedOfLight/sigma)*length/sqrt(z):0;
};
