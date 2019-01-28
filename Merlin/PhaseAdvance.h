/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef PhaseAdvance_h
#define PhaseAdvance_h 1

#include <string>
#include "AcceleratorModel.h"
#include "PSTypes.h"
#include "LatticeFunctions.h"
#include "TLAS.h"

using namespace TLAS;

/**
 * Class to calculate the phase advance of a given element
 * or the phase advance between two elements
 * or the transfer matrix between two elements
 * We assume that the lattice starts at the beginning of the AcceleratorModel
 */
class PhaseAdvance
{
public:
	PhaseAdvance(AcceleratorModel* aModel, LatticeFunctionTable* aTwiss, double refMomentum);

	void SetDelta(double new_delta);
	void ScaleBendPathLength(double scale);

	/**
	 * PA between two lattice elements
	 * can take either the ID number of the elements in the AcceleratorModel
	 * or the names of the two elements
	 */
	double PhaseAdvanceBetween(int n1, int n2, bool horizontal);
	double PhaseAdvanceBetween(string name1, string name2, bool horizontal);

	/**
	 * PA at a given element
	 * can take either the ID number of the element in the AcceleratorModel
	 * or the name of the element
	 */
	double PhaseAdvanceBetween(int n, bool horizontal);
	double PhaseAdvanceBetween(string name, bool horizontal);

	/**
	 * Calculates the transfer matrix between two lattice elements
	 */
	RealMatrix TransferMapBetween(int n1, int n2);

	// Simple functions for the user
	double GetPhaseAdvanceX(int n2, int n1 = 0);
	double GetPhaseAdvanceY(int n2, int n1 = 0);

	pair<double, double> CalcIntegerPart(int n);

private:
	AcceleratorModel* theModel;
	LatticeFunctionTable* theTwiss;
	double p0;
	double delta;
	double bendscale;
};

#endif
