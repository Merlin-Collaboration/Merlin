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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef Dispersion_h
#define Dispersion_h 1

#include <fstream>
#include "AcceleratorModel/AcceleratorModel.h"

class Dispersion
{
public:
    Dispersion(AcceleratorModel* aModel, double refMomentum);
    void FindDispersion(int n=0);
    void FindRMSDispersion(ofstream* file=0);
    double Dx, Dxp, Dy, Dyp;
    double DxRMS, DyRMS;
    double SetDelta(double new_delta);

private:
    AcceleratorModel* theModel;
    double p0;
    double delta;
};

#endif
