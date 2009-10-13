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

#ifndef TransferMatrix_h
#define TransferMatrix_h 1

#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/PSTypes.h"
#include "TLAS/TLAS.h"

using namespace TLAS;

class TransferMatrix
{
public:
    TransferMatrix(AcceleratorModel* aModel, double refMomentum);
    void FindTM(RealMatrix& M, PSvector& orbit);
    void FindClosedOrbitTM(RealMatrix& M, PSvector& orbit);
    void FindTM(RealMatrix& M);
    void Radiation(bool flag);
    void SetRadStepSize(double rad_stepsize);
    void SetRadNumSteps(int rad_numsteps);
    void ScaleBendPathLength(double scale);

    void SetObservationPoint(int n);
    void SetDelta(double new_delta);

private:
    AcceleratorModel* theModel;
    double p0;

    bool radiation;
    int obspnt;

    double delta;
    double radstepsize;
    int radnumsteps;
    double bendscale;

};

#endif
