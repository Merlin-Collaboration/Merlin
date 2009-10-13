/////////////////////////////////////////////////////////////////////////
// Class AcceleratorWithErrors
// Represents the physical accelerator with errors.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/13 11:46:56 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_AcceleratorWithErrors
#define _h_AcceleratorWithErrors

#include "Accelerator.h"
#include <list>

class AlignmentError;

class AcceleratorWithErrors : public Accelerator {
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
