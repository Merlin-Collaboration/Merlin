/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef TransferMatrix_h
#define TransferMatrix_h 1

#include "AcceleratorModel.h"
#include "PSTypes.h"
#include "TLAS.h"

using namespace TLAS;

class TransferMatrix
{
public:
	TransferMatrix(AcceleratorModel* aModel, double refMomentum);

	void FindTM(RealMatrix& M);
	void FindTM(RealMatrix& M, PSvector& orbit);
	void FindTM(RealMatrix& M, PSvector& orbit, int n1, int n2);

	void FindClosedOrbitTM(RealMatrix& M, PSvector& orbit);
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
