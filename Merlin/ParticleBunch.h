/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ParticleBunch_h
#define ParticleBunch_h 1

#include "merlin_config.h"
#include "PSTypes.h"
#include "Bunch.h"
#include "PhysicalConstants.h"

class Aperture; //#include "Aperture.h"
class ParticleDistributionGenerator; //#include "ParticleDistributionGenerator.h"
class BeamData; //#include "BeamData.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include <time.h>
#endif

class Transform3D;
using namespace PhysicalConstants;

namespace ParticleTracking
{
class ParticleBunchFilter; //#include "BunchFilter.h"
/**
 *	Representation of a particle.
 */
typedef PSvector Particle;

/**
 *	A Bunch which is represented by an ensemble of
 *	(macro-)particles.
 *
 *	The first particle in the bunch is the centroid.
 */
class ParticleBunch: public Bunch
{
public:

	typedef Particle particle_type;
	typedef PSvectorArray::iterator iterator;
	typedef PSvectorArray::const_iterator const_iterator;

	/**
	 *	Constructs a ParticleBunch using the specified momentum,
	 *	total charge and the particle array. Note that on exit,
	 *	particles is empty.
	 */
	ParticleBunch(double P0, double Q, PSvectorArray& particles);

	/**
	 *	Read phase space vectors from specified input stream.
	 */
	ParticleBunch(double P0, double Q, std::istream& is);

	/**
	 *	Constructs an empty ParticleBunch with the specified
	 *	momentum P0 and charge per macro particle Qm (default =
	 *	+1).
	 */
	ParticleBunch(double P0, double Qm = 1);

	/**
	 * Constructs an ParticleBunch with coordinates generated from a
	 * random distribution matched to a beam. Particles can be filtered
	 * using an optional ParticleBunchFilter.
	 */

	ParticleBunch(size_t np, const ParticleDistributionGenerator & generator, const BeamData& beam,
		ParticleBunchFilter* filter = nullptr);

	/**
	 *	Returns the total charge (in units of e).
	 *	@return Total charge (e)
	 */
	virtual double GetTotalCharge() const;

	virtual PSmoments& GetMoments(PSmoments& sigma) const;
	virtual PSmoments2D& GetProjectedMoments(PScoord u, PScoord v, PSmoments2D& sigma) const;
	virtual PSvector& GetCentroid(PSvector& p) const;
	virtual Point2D GetProjectedCentroid(PScoord u, PScoord v) const;

	/**
	 * Calculate mean (first) and rms (second) of the specified coordinate.
	 */
	std::pair<double, double> GetMoments(PScoord u) const;

	/**
	 *	Set the reference momentum to the mean (centroid)
	 *	momentum of the bunch. Returns the new value in GeV/c.
	 *
	 *	@return New reference momentum (GeV/c)
	 */
	virtual double AdjustRefMomentumToMean();

	/**
	 *	Adjusts the reference moment by the relative dpp.
	 */
	virtual double AdjustRefMomentum(double dpp);

	/**
	 *	Set the reference time to the mean (centroid) time of
	 *	the bunch. Returns the new value in ct (meters).
	 *
	 *	@return New reference time in ct (m)
	 */
	virtual double AdjustRefTimeToMean();

	/**
	 *	Used to generate a 1-D profile of the bunch projected
	 *	onto the specified coordinate. The total area of the
	 *	historgram is normalised to unity.
	 */
	virtual Histogram& ProjectDistribution(PScoord axis, Histogram& hist) const;

	/**
	 *	Apply the specified 3D coordinate transformation to the
	 *	bunch. Returns true if successful (note that it may not
	 *	be possible to apply a general 3D transformation to some
	 *	concrete bunch representations).
	 */
	virtual bool ApplyTransformation(const Transform3D& t);

	/**
	 *	Sorts the particles into ascending ct order.
	 */
	virtual void SortByCT();

	/**
	 *	Output a bunch-model dependent representation to the
	 *	specified stream.
	 */
	virtual void Output(std::ostream& os) const;
	virtual void Output(std::ostream& os, bool show_header) const;
	virtual void OutputIndexParticle(std::ostream& os, int index) const;
	virtual void Input(double Q, std::istream& is);

	/**
	 *	Add a (macro-)particle to the bunch.
	 */
	virtual size_t AddParticle(const Particle& p);

	/**
	 *	Sets the particle charge per macro-particle.
	 */
	void SetMacroParticleCharge(double q);

	ParticleBunch::iterator begin();
	ParticleBunch::iterator end();
	ParticleBunch::const_iterator begin() const;
	ParticleBunch::const_iterator end() const;
	size_t size() const;
	virtual void push_back(const Particle& p);
	virtual ParticleBunch::iterator erase(ParticleBunch::iterator p);
	void reserve(const size_t n);

	PSvectorArray& GetParticles();
	const PSvectorArray& GetParticles() const;

	/**
	 *	Returns the first particle in the bunch.
	 *	@return First particle in the bunch
	 */
	const Particle& FirstParticle() const;
	Particle& FirstParticle();

	/**
	 * Sets the centroid of the particle bunch to be equal to the zeroth
	 * particle, by moving each coordinate by the required amount.
	 */
	void SetCentroid();

	/**
	 * Sets the centroid of the particle bunch to be exactly x0 and updates
	 * the zeroth particle to x0.
	 */
	void SetCentroid(const Particle& x0);

	/**
	 * Removes all particles from the bunch
	 */
	void clear();

	/**
	 * Swaps particles with another ParticleBunch
	 */
	void swap(ParticleBunch newbunch);

	/**
	 * Init flag
	 */
	bool init;

	/**
	 * Number of coordinates involved in our Particle type
	 */
	int coords;

	/**
	 * Per-particle type scattering.
	 * This virtual function will be replaced in specific particle bunch classes
	 * to enable the relevant physics for each particle type (protons, electrons,
	 * etc).
	 */
	virtual int Scatter(Particle&, double length, const Aperture*)
	{
		return 0;
	}

	void SetScatterConfigured(bool);
	bool ScatterConfigured;
	size_t ScatteringPhysicsModel;
	double int_s;

	void SetIntS(double);

	double GetIntS()
	{
		return int_s;
	}

	/**
	 * Checks if the particle type is stable or not, returns true if the particle is considered stable.
	 *
	 * @retval true If particle is considered stable
	 * @retval false If particle is considered unstable
	 */
	virtual bool IsStable() const;

	/**
	 * Access method: Get particle mass
	 */
	virtual double GetParticleMass() const;

	/**
	 * Access method: Get particle mass (MeV)
	 */
	virtual double GetParticleMassMeV() const;

	/**
	 * Access method: Get particle lifetime
	 */
	virtual double GetParticleLifetime() const;

#ifdef ENABLE_MPI
	/**
	 * Destructor - cleans up MPI code
	 */
	~ParticleBunch();

	/**
	 * For send/recv buffers
	 */
	double* particle_send_buffer;
	double* particle_recv_buffer;

	/**
	 * For nanosecond timers
	 * timespec t_initial,t_final;
	 * double time_per_particle,t_delta;
	 */

	/**
	 * State information
	 */
	int MPI_size, MPI_rank;

	/**
	 * A particle
	 */
	MPI::Datatype MPI_Particle;

	/**
	 * MPI status
	 */
	MPI::Status MPI_status;

	/**
	 * Create particle type
	 */
	virtual void Create_MPI_particle();

	/**
	 * Init
	 */
	void MPI_Initialize();

	/**
	 * Gather to master
	 */
	void gather();

	/**
	 * push to nodes
	 */
	void distribute();

	/**
	 * node2master
	 */
	void master_recv_particles_from_nodes();
	/**
	 * master2node
	 */
	void node_recv_particles_from_master();

	void master_send_particles_to_nodes();
	void node_send_particles_to_master();

	/**
	 * Update reference momentum on master/nodes
	 */
	void SendReferenceMomentum();

	/**
	 * Finalize
	 */
	void MPI_Finalize();

	/**
	 * Check if MPI is active on this machine
	 */
	void Check_MPI_init();

	// For processor names
//	char proc_name[MPI_MAX_PROCESSOR_NAME];
//	int char_array_len;

#endif
private:

	/**
	 *	Charge per macro-particle
	 */
	double qPerMP;

protected:

	PSvectorArray pArray;

};
inline void ParticleBunch::swap(ParticleBunch newbunch)
{
	//cout << "Before " << size() << "\t" << newbunch.size() << endl;
	pArray.swap(newbunch.pArray);
	//cout << "After " << size() << "\t" << newbunch.size() << endl;
}

inline size_t ParticleBunch::AddParticle(const Particle& p)
{
	pArray.push_back(p);
	return size();
}

inline void ParticleBunch::SetMacroParticleCharge(double q)
{
	qPerMP = q;
	SetChargeSign(q);
}

inline ParticleBunch::iterator ParticleBunch::begin()
{
	return pArray.begin();
}

inline ParticleBunch::iterator ParticleBunch::end()
{
	return pArray.end();
}

inline void ParticleBunch::push_back(const Particle& p)
{
	AddParticle(p);
}

inline ParticleBunch::const_iterator ParticleBunch::begin() const
{
	return pArray.begin();
}

inline ParticleBunch::const_iterator ParticleBunch::end() const
{
	return pArray.end();
}

inline size_t ParticleBunch::size() const
{
	return pArray.size();
}

inline void ParticleBunch::reserve(const size_t n)
{
	pArray.reserve(n);
}

inline ParticleBunch::iterator ParticleBunch::erase(ParticleBunch::iterator p)
{
	return pArray.erase(p);
}

inline PSvectorArray& ParticleBunch::GetParticles()
{
	return pArray;
}

inline const PSvectorArray& ParticleBunch::GetParticles() const
{
	return pArray;
}

inline const Particle& ParticleBunch::FirstParticle() const
{
	return pArray.front();
}

inline Particle& ParticleBunch::FirstParticle()
{
	return pArray.front();
}

inline void ParticleBunch::clear()
{
	pArray.clear();
}

inline void ParticleBunch::SetScatterConfigured(bool state)
{
	ScatterConfigured = state;
}

inline void ParticleBunch::SetIntS(double step)
{
	int_s = step;
}

} // end namespace ParticleTracking

#endif
