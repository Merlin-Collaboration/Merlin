/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef PSvector_h
#define PSvector_h 1

#include "merlin_config.h"

#include <cstring>
#include <iostream>
#include <vector>
#include <cassert>
#include <cstdint>
#include "LinearAlgebra.h"

// Phase space coordinate indices
typedef int PScoord;
#define ps_X  0
#define ps_XP 1
#define ps_Y  2
#define ps_YP 3
#define ps_CT 4
#define ps_DP 5

#define PS_LENGTH 6

class PSvector
{
public:

	PSvector()
	{
	}

	explicit PSvector(double x)
	{
		std::fill(v, v + PS_LENGTH, x);
	}

	//	Component accessors.
	double x() const
	{
		return v[0];
	}
	double y() const
	{
		return v[2];
	}
	double ct() const
	{
		return v[4];
	}
	double xp() const
	{
		return v[1];
	}
	double yp() const
	{
		return v[3];
	}
	double dp() const
	{
		return v[5];
	}
	int64_t id() const
	{
		return pid;
	}

	//	Array access.
	double operator [](PScoord coord) const
	{
		assert(coord >= 0 && coord < PS_LENGTH);
		return v[coord];
	}

	//	Component mutators.
	double& x()
	{
		return v[0];
	}
	double& y()
	{
		return v[2];
	}
	double& ct()
	{
		return v[4];
	}
	double& xp()
	{
		return v[1];
	}
	double& yp()
	{
		return v[3];
	}
	double& dp()
	{
		return v[5];
	}
	int64_t& id()
	{
		return pid;
	}

	/**
	 *	Array access.
	 */
	double& operator [](PScoord coord)
	{
		assert(coord >= 0 && coord < PS_LENGTH);
		return v[coord];
	}

	/**
	 *	Conversion to a RealVector.
	 */
	operator RealVector() const
	{
		return RealVector(v, PS_LENGTH);
	}

	bool operator ==(const PSvector& psv) const
	{
		return memcmp(v, psv.v, PS_LENGTH * sizeof(double)) == 0;
	}

	bool operator !=(const PSvector& psv) const
	{
		return memcmp(v, psv.v, PS_LENGTH * sizeof(double)) != 0;
	}

	/**
	 *	Sets the vector to zero.
	 */
	void zero()
	{
		std::fill(v, v + PS_LENGTH, 0.0);
	}

	/**
	 *	Arithmetic assignment
	 */
	PSvector& operator +=(const PSvector& p)
	{
		double *q = v;
		const double *r = p.v;
		while(q != (v + PS_LENGTH))
		{
			*(q++) += *(r++);
		}
		return *this;
	}

	PSvector& operator -=(const PSvector& p)
	{
		double *q = v;
		const double *r = p.v;
		while(q != (v + PS_LENGTH))
		{
			*(q++) -= *(r++);
		}
		return *this;
	}

	PSvector& operator *=(double x)
	{
		for(double *q = v; q != v + PS_LENGTH; q++)
		{
			(*q) *= x;
		}
		return *this;
	}

	PSvector& operator /=(double x)
	{
		for(double *q = v; q != v + PS_LENGTH; q++)
		{
			(*q) /= x;
		}
		return *this;
	}

	/**
	 * binary arithmetic operators
	 */
	PSvector operator+(const PSvector& rhs) const
	{
		PSvector rv(*this);
		return rv += rhs;
	}

	PSvector operator-(const PSvector& rhs) const
	{
		PSvector rv(*this);
		return rv -= rhs;
	}

	/**
	 *	Reads the vector as six floating point numbers,
	 *	separated by spaces and terminated by a newline.
	 */
	friend std::ostream& operator<<(std::ostream& os, const PSvector& v);

	/**
	 *	Outputs the vector in row form, space delimited with a
	 *	terminating newline.
	 */
	friend std::istream& operator>>(std::istream& is, PSvector& v);

private:
	double v[PS_LENGTH];
	int64_t pid;
};

/**
 *	A linear array of PSvector objects.
 */
typedef std::vector<PSvector> PSvectorArray;

/**
 *	Representation of a particle.
 */
namespace ParticleTracking
{
typedef PSvector Particle;
}
#endif
