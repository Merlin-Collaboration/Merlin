/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ACG.h"
#include "Normal.h"
#include "Uniform.h"
#include "Poisson.h"
#include "Landau.h"
#include <cassert>
#include "RandomNG.h"

namespace
{

// table size for random number generator
#define TABLE_SIZE 100

}

RandGenerator* RandomNG::generator;

RandGenerator::RandGenerator (unsigned iseed)
	: nseed(iseed),
	  gen(nullptr),gaussGen(nullptr),uniformGen(nullptr),poissonGen(nullptr),landauGen(nullptr)
{
	reset(nseed);
}

RandGenerator::~RandGenerator ()
{
	if(gen)
	{
		delete gen;
	}
	if(gaussGen)
	{
		delete gaussGen;
	}
	if(uniformGen)
	{
		delete uniformGen;
	}
	if(poissonGen)
	{
		delete poissonGen;
	}
	if(landauGen)
	{
		delete landauGen;
	}
}



void RandGenerator::reset ()
{
	assert(gen);
	gen->reset();
	ResetGenerators();
}

void RandGenerator::reset (unsigned iseed)
{
	if(gen)
	{
		delete gen;
	}
	nseed = iseed;
	gen = new ACG(nseed,TABLE_SIZE);
	ResetGenerators();
}

double RandGenerator::normal (double mean, double variance)
{
	assert(gen && variance>=0);
	gaussGen->mean(mean);
	gaussGen->variance(variance);
	return (*gaussGen)();
}

double RandGenerator::normal (double mean, double variance, double cutoff)
{
	assert(gen);
	if(cutoff==0)
	{
		return normal(mean,variance);
	}

	cutoff=fabs(cutoff)*sqrt(variance);

	gaussGen->mean(mean);
	gaussGen->variance(variance);
	double x=(*gaussGen)();
	while(fabs(x-mean)>cutoff)
	{
		x=(*gaussGen)();
	}
	return x;
}

double RandGenerator::uniform (double low, double high)
{
	assert(gen);
	uniformGen->low(low);
	uniformGen->high(high);
	return (*uniformGen)();
}

double RandGenerator::poisson (double u)
{
	assert(gen);
	poissonGen->mean(u);
	return (*poissonGen)();
}

double RandGenerator::landau ()
{
	assert(gen);
	return (*landauGen)();
}

void RandGenerator::init (unsigned iseed)
{
	reset(iseed);
}

void RandGenerator::ResetGenerators ()
{
	if(gaussGen)
	{
		delete gaussGen;
	}
	if(uniformGen)
	{
		delete uniformGen;
	}
	if(poissonGen)
	{
		delete poissonGen;
	}
	if(landauGen)
	{
		delete landauGen;
	}

	gaussGen = new Normal(0,1,gen);
	uniformGen = new Uniform(0,1,gen);
	poissonGen = new Poisson(1,gen);
	landauGen = new Landau(gen);
}
