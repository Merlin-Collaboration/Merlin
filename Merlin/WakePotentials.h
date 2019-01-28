/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef WakePotentials_h
#define WakePotentials_h 1

#include "merlin_config.h"
#include "BunchProcess.h"

/**
 * Abstract class for calculating the longitudinal and
 * transverse single-bunch wakefield potentials (Greens
 * functions).
 */
class WakePotentials
{

public:

	WakePotentials(double r, double s) :
		csr(false), expectedProcess(nullptr), radius(r), conductivity(s)
	{
	}
	WakePotentials() :
		csr(false), expectedProcess(nullptr)
	{
	}                                                            // back to the original constructor
	//WakePotentials() : csr(false) {}   // back to the original constructor

	virtual ~WakePotentials()
	{
	}

	virtual double Wlong(double z) const = 0;
	virtual double Wtrans(double z) const = 0;

	bool Is_CSR() const
	{
		return csr;
	}

	BunchProcess* GetExpectedProcess() const
	{
		return expectedProcess;
	}

	void SetExpectedProcess(BunchProcess* p)
	{
		expectedProcess = p;
	}

protected:
	bool csr;

private:
	WakePotentials(const WakePotentials& wake);
	WakePotentials& operator=(const WakePotentials& wake);
	BunchProcess* expectedProcess;
	double radius;
	double conductivity;
};

#endif
