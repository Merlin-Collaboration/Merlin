/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorErrors_h
#define AcceleratorErrors_h 1

#include "AcceleratorModel.h"
#include "LatticeFrame.h"
#include "ComponentFrame.h"
#include "SequenceFrame.h"
#include "SupportStructure.h"

#include "StringPattern.h"
#include "NumericalConstants.h"
#include "RandomNG.h"

#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

/**
 * A class to apply random normally distributed errors to a beamline.
 * Both shift (transverse and longitudinal) and rotational errors can be generated.
 * @author Dirk Kruecker
 * @date 2008-12-01
 */
class AcceleratorErrors
{
public:

	/**
	 * Constructor
	 */
	AcceleratorErrors() :
		vx(0), vy(0), vz(0), mx(0), my(0), mz(0), clear(0), log(nullptr)
	{
	}

	/**
	 * Set the errors to apply to a Beamline.
	 * This function clears out any previous LatticeFrame movements that have been applied.
	 * @see AddErrors
	 * @param[in] xrms The root mean square value for x errors.
	 * @param[in] yrms The root mean square value for y errors.
	 * @param[in] zrms The root mean square value for z errors.
	 * @param[in] meanx The root mean square value for x errors.
	 * @param[in] meany The root mean square value for y errors.
	 * @param[in] meanz The root mean square value for z errors.
	 */
	void SetErrors(double xrms = 0, double yrms = 0, double zrms = 0, double meanx = 0, double meany = 0, double meanz =
		0)
	{
		vx = xrms * xrms;
		vy = yrms * yrms;
		vz = zrms * zrms;
		mx = meanx;
		my = meany;
		mz = meanz;
		clear = true;
	}

	/**
	 * Set the errors to apply to a Beamline.
	 * This function does NOT clear out any previous LatticeFrame movements that have been applied.
	 * It will add to any movements that have been applied.
	 * @see SetErrors
	 * @param[in] xrms The root mean square value for x errors.
	 * @param[in] yrms The root mean square value for y errors.
	 * @param[in] zrms The root mean square value for z errors.
	 * @param[in] meanx The root mean square value for x errors.
	 * @param[in] meany The root mean square value for y errors.
	 * @param[in] meanz The root mean square value for z errors.
	 */
	void AddErrors(double xrms = 0, double yrms = 0, double zrms = 0, double meanx = 0, double meany = 0, double meanz =
		0)
	{
		vx = xrms * xrms;
		vy = yrms * yrms;
		vz = zrms * zrms;
		mx = meanx;
		my = meany;
		mz = meanz;
		clear = false;
	}

	/**
	 * Apply the set shift errors to all elements matching a string pattern to a selected (sub) Beamline.
	 * @param[in] b The Beamline to apply the shift errors to.
	 * @param[in] p The string pattern for the names of elements to match for the application of errors.
	 */
	void ApplyShifts(AcceleratorModel::Beamline& b, const string& p);

	/**
	 * Apply the set rotational errors to all elements matching a string pattern to a selected (sub) Beamline.
	 * @param[in] b The Beamline to apply the rotational errors to.
	 * @param[in] p The string pattern for the names of elements to match for the application of errors.
	 */
	void ApplyRotations(AcceleratorModel::Beamline& b, const string& p);

	/**
	 * Sets the log stream to output logging information to.
	 * @param[in] l The stream to use for output.
	 */
	void Report(ostream* l)
	{
		log = l;
	}

private:

	double vx, vy, vz;
	double mx, my, mz;
	bool clear;
	ostream* log;

	//Copy protection
	AcceleratorErrors(const AcceleratorErrors& rhs);
	AcceleratorErrors& operator=(const AcceleratorErrors& rhs);
};
#endif
