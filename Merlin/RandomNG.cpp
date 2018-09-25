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
#include "LandauDistribution.h"

namespace
{

// table size for random number generator
#define TABLE_SIZE 100

}

RandGenerator::RandGenerator(unsigned iseed) :
	nseed(iseed), gen(nullptr), gaussGen(nullptr), uniformGen(nullptr), poissonGen(nullptr), landauGen(nullptr)
{
	reset(nseed);
}

RandGenerator::~RandGenerator()
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

void RandGenerator::reset()
{
	assert(gen);
	gen->reset();
	ResetGenerators();
}

void RandGenerator::reset(unsigned iseed)
{
	if(gen)
	{
		delete gen;
	}
	nseed = iseed;
	gen = new ACG(nseed, TABLE_SIZE);
	ResetGenerators();
}

double RandGenerator::normal(double mean, double variance)
{
	assert(gen && variance >= 0);
	gaussGen->mean(mean);
	gaussGen->variance(variance);
	return (*gaussGen)();
}

double RandGenerator::normal(double mean, double variance, double cutoff)
{
	assert(gen);
	if(cutoff == 0)
	{
		return normal(mean, variance);
	}

	cutoff = fabs(cutoff) * sqrt(variance);

	gaussGen->mean(mean);
	gaussGen->variance(variance);
	double x = (*gaussGen)();
	while(fabs(x - mean) > cutoff)
	{
		x = (*gaussGen)();
	}
	return x;
}

double RandGenerator::uniform(double low, double high)
{
	assert(gen);
	uniformGen->low(low);
	uniformGen->high(high);
	return (*uniformGen)();
}

double RandGenerator::poisson(double u)
{
	assert(gen);
	poissonGen->mean(u);
	return (*poissonGen)();
}

double RandGenerator::landau()
{
	assert(gen);
	return (*landauGen)();
}

void RandGenerator::init(unsigned iseed)
{
	reset(iseed);
}

void RandGenerator::ResetGenerators()
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

	gaussGen = new Normal(0, 1, gen);
	uniformGen = new Uniform(0, 1, gen);
	poissonGen = new Poisson(1, gen);
	landauGen = new Landau(gen);
}

std::vector<std::uint32_t> RandomNG::master_seed;
std::unique_ptr<std::mt19937_64> RandomNG::generator;

std::unordered_map<size_t, std::mt19937_64> RandomNG::generator_store;

void RandomNG::init()
{
	const int nseeds = 8;
	master_seed.clear();
	master_seed.reserve(nseeds);
	for(int i = 0; i < nseeds; i++)
	{
		master_seed.push_back(std::random_device{} ());
	}
	reset();
}

void RandomNG::init(std::uint32_t iseed)
{
	master_seed = {iseed};
	reset();
}

void RandomNG::init(const std::vector<std::uint32_t>& iseed)
{
	master_seed = iseed;
	reset();
}

void RandomNG::reset()
{
	std::seed_seq ss(master_seed.begin(), master_seed.end());
	generator.reset(new std::mt19937_64{ss});
}

void RandomNG::reset(std::uint32_t iseed)
{
	master_seed = {iseed};
	reset();
}

void RandomNG::reset(const std::vector<std::uint32_t>& iseed)
{
	master_seed = iseed;
	reset();
}

//std::vector<std::uint32_t> RandomNG::getSeed()
std::uint32_t RandomNG::getSeed()
{
	return master_seed[0];
}

double RandomNG::normal(double mean, double variance)
{
	if(!generator)
		not_seeded();
	std::normal_distribution<double> dist{mean, sqrt(variance)};
	return dist(*generator);
}
double RandomNG::normal(double mean, double variance, double cutoff)
{
	if(!generator)
	{
		not_seeded();
	}
	if(cutoff == 0)
	{
		return normal(mean, variance);
	}

	cutoff = fabs(cutoff) * sqrt(variance);
	std::normal_distribution<double> dist{mean, sqrt(variance)};

	double x;
	do
	{
		x = dist(*generator);
	} while(fabs(x - mean) > cutoff);
	return x;
}
double RandomNG::uniform(double low, double high)
{
	if(!generator)
	{
		not_seeded();
	}
	std::uniform_real_distribution<double> dist{low, high};
	return dist(*generator);
}
double RandomNG::poisson(double u)
{
	if(!generator)
	{
		not_seeded();
	}
	std::poisson_distribution<int> dist{u};
	return dist(*generator);
}
double RandomNG::landau()
{
	if(!generator)
	{
		not_seeded();
	}
	landau_distribution<double> dist{};
	return dist(*generator);
}

std::mt19937_64& RandomNG::getGenerator()
{
	return *generator;
}

std::mt19937_64& RandomNG::getLocalGenerator(size_t name_hash)
{
	auto search = generator_store.find(name_hash);
	if(search != generator_store.end())
	{
		return search->second;
	}
	else
	{
		// extend master seed
		std::vector<std::uint32_t> new_seed{master_seed};
		new_seed.push_back(name_hash);
		std::seed_seq ss(new_seed.begin(), new_seed.end());
		// create and store new generator
		auto new_gen = std::mt19937_64{ss};
		generator_store[name_hash] = new_gen;
		return generator_store[name_hash];
	}
}

void RandomNG::resetLocalGenerator(size_t name_hash)
{
	std::vector<std::uint32_t> new_seed{master_seed};
	new_seed.push_back(name_hash);
	std::seed_seq ss(new_seed.begin(), new_seed.end());
	// create and store new generator
	auto new_gen = std::mt19937_64{ss};
	generator_store[name_hash] = new_gen;
}

std::uint32_t hash_string(std::string s)
{
	return std::hash<std::string>{} (s);
}
