/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Measurement_h
#define Measurement_h 1

/**
 *	POD representing a physically measured quantity, which
 *	has a value and an associated error.
 */
struct Measurement
{
	Measurement(double v, double err);
	Measurement();

	double value;
	double error;

};

/**
 * Class Measurement
 */
inline Measurement::Measurement(double v, double err) :
	value(v), error(err)
{
}

inline Measurement::Measurement() :
	value(0), error(0)
{
}

#endif
