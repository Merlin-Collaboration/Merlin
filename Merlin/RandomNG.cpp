/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "RandomNG.h"
#include "LandauDistribution.h"

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

const std::vector<std::uint32_t>& RandomNG::getSeed()
{
	return master_seed;
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
	if(!generator)
	{
		not_seeded();
	}
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
	if(!generator)
	{
		not_seeded();
	}
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
