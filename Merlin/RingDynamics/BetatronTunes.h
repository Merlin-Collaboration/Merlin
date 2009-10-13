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
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef BetatronTunes_h
#define BetatronTunes_h 1

#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/PSTypes.h"

class BetatronTunes
{
public:
    BetatronTunes(AcceleratorModel* aModel, double refMomentum);
    void FindTunes(PSvector& particle, int ntrack = 256, bool diffusion = true);
    double Qx, Qy, dQx, dQy;

private:
    AcceleratorModel* theModel;
    double p0;
    double FindTune(vector<double>& data);
    void FFT(vector<double>& data);
    double amp(double a, double b, double c);
};

#endif
