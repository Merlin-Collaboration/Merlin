/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef LatticeFunctions_h
#define LatticeFunctions_h 1

#include "PSvector.h"

class LatticeFunction
{
public:
	LatticeFunction(int _i, int _j, int _k);
	~LatticeFunction();
	void GetIndices(int& _i, int& _j, int& _k);
	void AppendValue(double v);
	void ClearValues();
	double GetValue(int n);
	void Derivative(LatticeFunction* lfnM, LatticeFunction* lfnP, double dp);
	vector<double>::iterator begin();
	vector<double>::iterator end();
	int size();

private:
	int i, j, k;
	vector<double> value;
};

typedef vector<LatticeFunction*> vectorlfn;

class LatticeFunctionTable
{
public:
	LatticeFunctionTable(AcceleratorModel* aModel, double refMomentum);
	~LatticeFunctionTable();
	void AddFunction(int i, int j, int k);
	void UseDefaultFunctions();
	void UseOrbitFunctions();
	void RemoveFunction(int i, int j, int k);
	void RemoveAllFunctions();
	void Calculate(PSvector* p = nullptr, RealMatrix* M = nullptr);
	void CalculateEnergyDerivative();
	double Value(int i, int j, int k, int ncpt);
	void PrintTable(ostream& os, int n1 = 0, int n2 = -1);
	void Size(int& rows, int& cols);
	int GetSPosIndex(double s);
	void SetDelta(double new_delta);
	void MakeTMSymplectic(bool flag);
	int NumberOfRows();
	void ScaleBendPathLength(double scale);
	double Mean(int i, int j, int k, int n1 = 0, int n2 = -1);
	double RMS(int i, int j, int k, int n1 = 0, int n2 = -1);

private:
	AcceleratorModel* theModel;
	double p0;
	double delta;
	double bendscale;
	bool symplectify;
	bool orbitonly;

	vectorlfn lfnlist;

	double DoCalculate(double cscale = 0, PSvector* pInit = nullptr, RealMatrix* MInit = nullptr);
	double DoCalculateOrbitOnly(double cscale = 0, PSvector* pInit = nullptr);
	vectorlfn::iterator GetColumn(int i, int j, int k);
};

#endif
