/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RandomNG_h
#define RandomNG_h 1

#include "merlin_config.h"
#include <random>
#include <memory>
#include <cstdint>
#include <iostream>
#include <unordered_map>

/**
 * Singleton class for generating continuous floating point numbers from specific distributions.
 * * normal (Gaussian)
 * * uniform distributions
 * * poisson
 * * landau
 *
 * Also provides access to the generator for more optimised usage.
 */

class RandomNG
{
public:
	RandomNG() = delete;

	/**
	 * Initialise the generator. One form of init() be called before any
	 * other generator function. With no argument a seed is automatically
	 * generated using the system random source.
	 */
	static void init();

	/// Initialise with single unsigned int as seed
	static void init(std::uint32_t iseed);
	/// Initialise with list of unsigned ints as seed
	static void init(const std::vector<std::uint32_t>& iseed);

	/**
	 * Resets the seed for the generators to the last supplied
	 * seed value.
	 */
	static void reset();

	/**
	 * Reset the generator with a new seed
	 */
	static void reset(std::uint32_t iseed);
	/// Reset the generator with a new seed
	static void reset(const std::vector<std::uint32_t>& iseed);

	/**
	 * Get the seed for the generator.
	 */
	static const std::vector<std::uint32_t>& getSeed();

	/**
	 * Generates a random number from a normal (Gaussian)
	 * distribution with the specified mean and variance.
	 */
	static double normal(double mean, double variance);

	/**
	 * Generates a random number from a normal (Gaussian)
	 * distribution with the specified mean and variance. The
	 * resulting distribution is truncated to +/-cutoff
	 * standard deviations.
	 */
	static double normal(double mean, double variance, double cutoff);

	/**
	 * Generates a uniform random number in the range
	 * |low,high> inclusive.
	 */
	static double uniform(double low, double high);

	/**
	 * Generates a Poisson random number in with a given expected value.
	 */
	static double poisson(double u);

	/**
	 * Generates a  random number.
	 */
	static double landau();

	/// Gives a reference to the actual generator
	static std::mt19937_64& getGenerator();

	/**
	 * Get a new generator to be used within a physics process class
	 *
	 * This can improve reproducibility, as each local generator is independent.
	 * name_hash is added to the list of seeds, and should be set using
	 * the name of the calling class, e.g.
	 *
	 *     auto gen = RandomNG::getLocalGenerator(hash_string("MyPhysicsProcess"));
	 *
	 * The generator is held within RandomNG::generator_store, so the same
	 * generator will be returned if called with the same name_hash
	 */
	static std::mt19937_64& getLocalGenerator(size_t name_hash);

	/// Reset a given local generator
	static void resetLocalGenerator(size_t name_hash);

private:
	static std::vector<std::uint32_t> master_seed;
	static std::unique_ptr<std::mt19937_64> generator;

	static std::unordered_map<size_t, std::mt19937_64> generator_store;

	static void not_seeded()
	{
		std::cerr << "WARN: Random number generator not initiated, using auto seeding" << std::endl;
		std::cerr << "WARN: It is recommended to call MerlinRandom::init() with a seed" << std::endl;
		init();
	}
};

std::uint32_t hash_string(std::string s);
#endif
