/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_AcceleratorWithErrors
#define _h_AcceleratorWithErrors

#include "Accelerator.h"
#include <list>

class AlignmentError;

class AcceleratorWithErrors: public Accelerator
{
public:

	AcceleratorWithErrors(const std::string& name, AcceleratorModel*, BeamData*);

	// RMS gaussian error definitions
	void TransverseErrors(const string& pattern, double x_rms, double y_rms, bool clear_transform);
	void RotationErrors(const string& pattern, double x_rms, double y_rms, double z_rms, bool clear_transform);
	void MagnetScaleError(const string& pattern, double rms);
	void KlystronErrors(const string& pattern, double v_rms, double phi_rms);
	void BPMresolution(double rms);
	void BPMlinearError(double rms);
	void InitialBeamJitterInSigma(double xrms, double yrms);

	void ApplyStaticErrors();

private:
	std::list<AlignmentError*> alignErrorDefs;
	double bpmCalErr2;
};
#endif
