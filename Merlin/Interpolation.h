/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Interpolation_h
#define Interpolation_h 1

#include "merlin_config.h"
#include "MerlinException.h"
#include "Range.h"
#include <vector>

/**
 * class Interpolation
 * An interpolation functor which interpolates values from a data table.
 * Currently only linear interpolation is assumed.
 */

class Interpolation
{
public:

	/**
	 * exception
	 */
	class BadRange: public MerlinException
	{
	public:
		BadRange(double x, const FloatRange& r);
		double value;
		FloatRange valid_range;
	};

	/**
	 * Implementation method for interpolation
	 */
	class Method
	{
	public:
		virtual ~Method()
		{
		}
		virtual double ValueAt(double x) const = 0;
	};

	/**
	 * Interpolation of equally spaced data points
	 */
	Interpolation(const std::vector<double>& yvals, double xmin, double dx);

	/**
	 * Interpolation of arbitrary spaced data points
	 */
	Interpolation(const std::vector<double>& xvals, const std::vector<double>& yvals);

	~Interpolation();

	double operator()(double x) const
	{
		return itsMethod->ValueAt(x);
	}

private:

	Method* itsMethod;
};

#endif
