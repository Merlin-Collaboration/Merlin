#ifndef ParticleDistributionGenerator_h
#define ParticleDistributionGenerator_h 1

#include "PSvector.h"
#include "RandomNG.h"

inline double RandomGauss(double variance, double cutoff)
{
	return cutoff==0 ? RandomNG::normal(0,variance) :  RandomNG::normal(0,variance,cutoff);
}

/**
* Base class for distribution generators. These can be used by
* ParticleTracking::ParticleBunch::ParticleBunch to construct bunches with
* a given distribution.
*
* Derived classes must override GenerateFromDistribution(), with a function
* that returns a single PSvector from the distribution.
*
*  Additional parameters can be passed to the constructors of derived classes.
*/
class ParticleDistributionGenerator
{
public:
	/**
	 * Returns a single PSvector from the distribution
	 */
	virtual PSvector GenerateFromDistribution() const = 0;
	virtual ~ParticleDistributionGenerator() {};
};

/**
 * Generator for a normal (Gaussian) distribution.
 */
class NormalParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	/**
	 * @param cutoffs_ Vector of cut off points in the distribution in each coordinate.
	 * Default zero gives no cut off.
	 */
	NormalParticleDistributionGenerator(PSvector cutoffs_ = PSvector(0)): cutoffs(cutoffs_) {};
	/**
	 * @param cutoff Cut off point in distribution, same in each coordinate
	 */
	NormalParticleDistributionGenerator(double cutoff): cutoffs(PSvector(cutoff)) {};
	virtual PSvector GenerateFromDistribution() const override;
private:
	PSvector cutoffs;
};

/**
 * Generator for a uniform (flat) distribution.
 */
class UniformParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	virtual PSvector GenerateFromDistribution() const override;
};

/**
 * Generator for a ring distribution.
 */
class RingParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	virtual PSvector GenerateFromDistribution() const override;
};

#endif
