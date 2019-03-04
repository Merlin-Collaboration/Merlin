/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_NormalTransform
#define _h_NormalTransform

#include "LinearAlgebra.h"
#include "BeamData.h"
#include "PSTypes.h"

// utility global functions for generating beam phase space
// normalisations

RealMatrix CouplingMatrix(double a, double b, double c, double d);
RealMatrix AlphaMatrix(double ax, double ay);
RealMatrix BetaMatrix(double bx, double by);
RealMatrix DispersionMatrix(double Dx, double Dxp, double Dy, double Dyp);
RealMatrix NormalTransform(const BeamData& t);
RealMatrix DecoupleSigma(SigmaMatrix& S);
RealMatrix InverseBetaTransform(double bx, double by, double ax, double ay);

PSmoments& BeamDataToSigmaMtrx(const BeamData& t, PSmoments& S);
inline PSmoments BeamDataToSigmaMtrx(const BeamData& t)
{
	PSmoments S;
	return BeamDataToSigmaMtrx(t, S);
}

BeamData& SigmaMatrixToBeamData(const PSmoments& S0, BeamData& t);

double ProjectedEmittance(const PSmoments& s, PScoord x1, PScoord x2);

pair<double, double> NormalModeEmittance(const PSmoments& S);

#endif
