/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ParticleBunchConstructor_h
#define ParticleBunchConstructor_h 1

#include "merlin_config.h"
#include "BeamData.h"
#include "BunchConstructor.h"
#include "MatrixMaps.h"
#include "ParticleBunch.h"
#include "ParticleBunchTypes.h"
#include "BunchFilter.h"
#include <typeinfo>
#include "PSvector.h"

namespace ParticleTracking
{

typedef enum
{
	normalDistribution,
	flatDistribution,
	ringDistribution,
	skewHaloDistribution,
	horizontalHaloDistribution1,
	verticalHaloDistribution1,
	horizontalHaloDistribution2,
	verticalHaloDistribution2

} DistributionType;

/**
 * Note this class is deprecated and should no longer be used. See \ref APIChanges .
 *
 * Constructs a particle bunch with random particles taken
 * from a 6D distribution. The phase space moments are
 * supplied as a BeamData struct. The form of the
 * distribution can be either normal (Gaussian), in which
 * case the data supplied in the BeamData struct is
 * interpreted as  RMS values, or flat, in which case the
 * data specifies the +- extents for a flat (rectangular)
 * distribution.
 *
 * An additional ParticleBunchFilter can be specified which
 * can further be used to modify the phase space
 * distribution generated.
 */
class ParticleBunchConstructor: public BunchConstructor
{
public:

	/**
	 * Constructor taking the beam data and the number of
	 * particles to generate.
	 */
	[[deprecated("Use the ParticleBunch constructor directly.")]] ParticleBunchConstructor(const BeamData& beam, size_t
		npart, DistributionType dist = normalDistribution);

	~ParticleBunchConstructor();

	/**
	 * Sets the bunch information that will be used to generate the particle bunch.
	 */
	void SetBunchData(const BeamData& beam);

	/**
	 * Sets the number of particles to be generated.
	 */
	void SetNumParticles(size_t npart);

	/**
	 * Sets the distribution cut-off in standard deviations.
	 * Passing 0 indicates that no cut-off be applied
	 * (default). Only applies to normal distributions. The
	 * second form allows individual cutoffs for each of the
	 * six phase-space coordinates.
	 */
	void SetDistributionCutoff(double cut);
	void SetDistributionCutoff(const PSvector& cut);

	/**
	 * Constructs a new ParticleBunch based on the current bunch parameter settings. Each call to ConstructBunch
	 * generates a new random distribution (seed). The bunch Index is ignored in this case.
	 * Required due to pure virtual in Bunch.h
	 */
	virtual Bunch* ConstructBunch(int bunchIndex = 0) const;

	//template <class T_bunch> T_bunch* ConstructParticleBunch () const;

	virtual void ConstructBunchDistribution(int bunchIndex = 0) const;
	PSvector GenerateFromDistribution() const;

	/**
	 * Sets the filter to be used during bunch construction. A NULL pointer indicates no filter.
	 */
	void SetFilter(ParticleBunchFilter* filter);

	//	Returns typed particle bunch.
	ParticleBunch* ConstructParticleBunch() const;
	template<class T_bunch> T_bunch* ConstructParticleBunch() const
	{
		ConstructBunchDistribution();
		return new T_bunch(beamdat.p0, beamdat.charge, pbunch);
	}

	/**
	 * If fc==true, causes the centroids (means) of the
	 * particle distribution to be set to the specified
	 * centroid. When fc==false (default), the actually
	 * generated centroids will contain statistical variations
	 * from the desired values.
	 */
	void ForceCentroid(bool fc);

private:
	void DoForceCentroid() const;

	size_t np;
	DistributionType dtype;
	PSvector cutoffs;
	BeamData beamdat;
	ParticleBunchFilter* itsFilter;

	/**
	 * The map used to transform a spherically symmetric
	 * distribution into the desired bunch distribution.
	 */
	RMtrx M;
	bool force_c;

	//Moved the pbunch to be a class member so that the bunch creation can be split up between distribution generation and "bunch" generation.
	mutable PSvectorArray pbunch;
};

inline void ParticleBunchConstructor::SetFilter(ParticleBunchFilter* filter)
{

	if(itsFilter)
	{
		delete itsFilter;
	}
	itsFilter = filter;

}

inline ParticleBunch* ParticleBunchConstructor::ConstructParticleBunch() const
{
	return ParticleBunchConstructor::ConstructParticleBunch<ParticleBunch>();
//    return static_cast<ParticleBunch*>(ParticleBunchConstructor::ConstructBunch<ParticleBunch>());
}

} //end namespace ParticleTracking

#endif
