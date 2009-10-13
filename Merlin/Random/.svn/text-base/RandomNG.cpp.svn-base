/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include "Random/ACG.h"
#include "Random/Normal.h"
#include "Random/Uniform.h"
#include "Random/Poisson.h"
#include <cassert>
#include "Random/RandomNG.h"

namespace {

// table size for random number generator
	#define TABLE_SIZE 100

};

RandGenerator* RandomNG::generator;

RandGenerator::RandGenerator (unsigned iseed)
        : nseed(iseed),
        gen(0),gaussGen(0),uniformGen(0),poissonGen(0)
{
    reset(nseed);
}

RandGenerator::~RandGenerator ()
{
    if(gen)
        delete gen;
    if(gaussGen)
        delete gaussGen;
    if(uniformGen)
        delete uniformGen;
    if(poissonGen)
        delete poissonGen;
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
        delete gen;
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
        return normal(mean,variance);

    cutoff=fabs(cutoff)*sqrt(variance);

    gaussGen->mean(mean);
    gaussGen->variance(variance);
    double x=(*gaussGen)();
    while(fabs(x-mean)>cutoff)
        x=(*gaussGen)();
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

void RandGenerator::init (unsigned iseed)
{
    reset(iseed);
}

void RandGenerator::ResetGenerators ()
{
    if(gaussGen)
        delete gaussGen;
    if(uniformGen)
        delete uniformGen;
    if(poissonGen)
        delete poissonGen;

    gaussGen = new Normal(0,1,gen);
    uniformGen = new Uniform(0,1,gen);
    poissonGen = new Poisson(1,gen);
}

