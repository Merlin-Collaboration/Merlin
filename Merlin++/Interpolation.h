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
#include <iostream>

/**
 * class Interpolation
 * An interpolation functor which interpolates values from a data table.
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

		const char* what() const noexcept
		{
			return msg.c_str();
		}
	private:
		std::string msg;
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
		virtual double ValueAt(double x,double *err=0x0) = 0;
		double evaluatepoly(double x,std::vector<double> &p){ 
			int n=p.size(); 
			double s=p[n-1];
			if(n>1) for(int j=2;j<=n;j++) s =s*x+p[n-j]; 
			return s;
		} 
	int itsOrder;
	};

	/**
	 * Interpolation of equally spaced data points
	 */
	Interpolation(const std::vector<double>& yvals, double xmin, double dx,int order=1);

	/**
	 * Interpolation of arbitrary spaced data points
	 */
	Interpolation(const std::vector<double>& xvals, const std::vector<double>& yvals,int order=1);
	Interpolation(const double* xvals, const double* yvals,int n,int order=1);

	~Interpolation();

	double operator()(double x) const
	{
		return itsMethod->ValueAt(x);
	}



	Method* itsMethod;
};

#endif
