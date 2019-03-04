/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>

#include "AcceleratorModel.h"

#include "ParticleBunch.h"
#include "ParticleTracker.h"
#include "RingDeltaTProcess.h"

#include "ClosedOrbit.h"
#include "TransferMatrix.h"
#include "PhaseAdvance.h"
#include "LatticeFunctions.h"

#include "MatrixPrinter.h"

using namespace ParticleTracking;

PhaseAdvance::PhaseAdvance(AcceleratorModel* aModel, LatticeFunctionTable* aTwiss, double refMomentum) :
	theModel(aModel), theTwiss(aTwiss), p0(refMomentum), delta(1.0E-8), bendscale(1E-16)
{
	//~ std::cout << "\n\tPhaseAdvance Class created" << std::endl;
	//~ std::cout << "please note that the following functions must be added to your LatticeFunctionTable: (0,0,1), (0,0,2), (0,0,3), using the following function:" << std::endl;
	//~ std::cout << "\n\t\t myLatticeFunctionTable->AddFunction(0,0,1);" << std::endl;
	//~ std::cout << "\n\t\t myLatticeFunctionTable->AddFunction(0,0,2);" << std::endl;
	//~ std::cout << "\n\t\t myLatticeFunctionTable->AddFunction(0,0,3);" << std::endl;
}

void PhaseAdvance::SetDelta(double new_delta)
{
	delta = new_delta;
}

void PhaseAdvance::ScaleBendPathLength(double scale)
{
	bendscale = scale;
}

double PhaseAdvance::PhaseAdvanceBetween(int n1, int n2, bool horizontal)
{

	double deltamu;

	pair<double, double> element1 = CalcIntegerPart(n1);
	pair<double, double> element2 = CalcIntegerPart(n2);

	if(horizontal)
	{
		deltamu = element2.first - element1.first;
	}
	else
	{
		deltamu = element2.second - element1.second;
	}

	//OBSOLETE/NOT WORKING
	//~ RealMatrix M = TransferMapBetween(n1,n2);

	//~ if(horizontal){
	//~ deltamu = asin( M(1,0) / (sqrt(theTwiss->Value(1,1,1,n1) * theTwiss->Value(1,1,1,n2))) );
	//~ }
	//~ else{
	//~ deltamu = asin( M(2,3) / (sqrt(theTwiss->Value(3,3,2,n1) * theTwiss->Value(3,3,2,n2))) );
	//~ }

	return deltamu;
}

double PhaseAdvance::PhaseAdvanceBetween(string name1, string name2, bool horizontal)
{
	int n1 = theModel->FindElementLatticePosition(name1);
	int n2 = theModel->AcceleratorModel::FindElementLatticePosition(name2);

	double deltamu = PhaseAdvance::PhaseAdvanceBetween(n1, n2, horizontal);

	return deltamu;
}

double PhaseAdvance::PhaseAdvanceBetween(int n, bool horizontal)
{
	double deltamu = PhaseAdvanceBetween(0, n, horizontal);

	return deltamu;
}
double PhaseAdvance::PhaseAdvanceBetween(string name, bool horizontal)
{
	int n = theModel->FindElementLatticePosition(name);
	double deltamu = PhaseAdvanceBetween(n, horizontal);

	return deltamu;
}

RealMatrix PhaseAdvance::TransferMapBetween(int n1, int n2)
{
	Particle p2(0.);

	if(n1 > n2)
	{
		std::cout
			<<
			"\n\tPhaseAdvance::TransferMapBetween: WARNING: n1>n2, tracker cannot compute \nt Use GetPhaseAdvanceX(n2, n1)"
			<< endl;
	}

	ClosedOrbit co(theModel, p0);
	co.SetDelta(delta);
	co.TransverseOnly(false);
	co.ScaleBendPathLength(bendscale);

	RealMatrix M(6);
	TransferMatrix tm(theModel, p0);
	tm.SetDelta(delta);
	tm.ScaleBendPathLength(bendscale);

	co.FindClosedOrbit(p2);
	//~ std::cout << "PhaseAdvance::TransferMapBetween: ClosedOrbit done" << endl;
	tm.FindTM(M, p2, n1, n2);
	//~ std::cout << "PhaseAdvance::TransferMapBetween: TransferMatrix done" << endl;

	return M;
}

double PhaseAdvance::GetPhaseAdvanceX(int n2, int n1)
{
	return PhaseAdvanceBetween(n1, n2, 1);
}

double PhaseAdvance::GetPhaseAdvanceY(int n2, int n1)
{
	return PhaseAdvanceBetween(n1, n2, 0);
}

pair<double, double> PhaseAdvance::CalcIntegerPart(int n)
{
	//Fractional Phase Advance stored in Twiss
	//MuX = theTwiss->Value(0,0,1,n)
	//MuY = theTwiss->Value(0,0,2,n)
	double mux(0.), muy(0.);
	int intmux(0), intmuy(0);
	double last_mux(0.), last_muy(0.);

	//Iterate through all elements to sum integer parts
	for(int i = 1; i < n; ++i)
	{
		if((theTwiss->Value(0, 0, 1, i) == last_mux) || (theTwiss->Value(0, 0, 1, i) <= 0.))
		{
		}
		else if((theTwiss->Value(0, 0, 1, i) < last_mux))
		{
			//~ if( (theTwiss->Value(0,0,1,i)/last_mux) < 0.9999999999){ //to stop iteration for small fluctuations
			if((theTwiss->Value(0, 0, 1, i) / last_mux) < 0.9999)  //to stop iteration for small fluctuations
			{
				++intmux;
			}
		}
		if(theTwiss->Value(0, 0, 2, i) == last_muy || (theTwiss->Value(0, 0, 2, i) <= 0.))
		{
		}
		else if(theTwiss->Value(0, 0, 2, i) < last_muy)
		{
			//~ if( (theTwiss->Value(0,0,2,i)/last_muy) < 0.9999999999){ //to stop iteration for small fluctuations
			if((theTwiss->Value(0, 0, 2, i) / last_muy) < 0.9999)  //to stop iteration for small fluctuations
			{
				++intmuy;
			}
		}
		last_mux = theTwiss->Value(0, 0, 1, i);
		last_muy = theTwiss->Value(0, 0, 2, i);
	}

	mux = intmux + theTwiss->Value(0, 0, 1, n);
	muy = intmuy + theTwiss->Value(0, 0, 2, n);
	//~ cout << "mux = " << mux << endl;
	//~ cout << "muy = " << muy << endl;

	return make_pair(mux, muy);
}
