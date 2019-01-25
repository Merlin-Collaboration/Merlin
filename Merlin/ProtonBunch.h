/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ProtonBunch_h
#define ProtonBunch_h 1

#include "ParticleBunch.h"
#include <iostream>
#include "Collimator.h"
#include "PhysicalConstants.h"
#include <vector>
#include "ElasticScatter.h"
#include "DiffractiveScatter.h"
#include "MerlinProfile.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
namespace ParticleTracking
{

class ProtonBunch: public ParticleBunch
{
	static const int ntally = 6;
	int tally[ntally];

public:

	/**
	 * Constructs a ProtonBunch using the specified momentum,
	 * total charge and the particle array. Note that on exit,
	 * particles is empty.
	 */
	ProtonBunch(double P0, double Q, PSvectorArray& particles) :
		ParticleBunch(P0, Q, particles), GotElastic(false), GotDiffractive(false)
	{
		SetUpProfiling();
	}

	/**
	 * Read phase space vectors from specified input stream.
	 */
	//ProtonBunch (double P0, double Q, std::istream& is) : ParticleBunch(P0, Q, is) {rng();};
	ProtonBunch(double P0, double Q, std::istream& is) :
		ParticleBunch(P0, Q, is), GotElastic(false), GotDiffractive(false)
	{
		SetUpProfiling();
	}

	/**
	 * Constructs an empty ProtonBunch with the specified
	 * momentum P0 and charge per macro particle Qm (default =
	 *  +1).
	 */
	//ProtonBunch (double P0, double Qm = 1) : ParticleBunch(P0, Qm) {rng();};
	ProtonBunch(double P0, double Qm = 1) :
		ParticleBunch(P0, Qm), GotElastic(false), GotDiffractive(false)
	{
		SetUpProfiling();
	}

	ProtonBunch(size_t np, const ParticleDistributionGenerator & generator, const BeamData& beam,
		ParticleBunchFilter* filter = nullptr) :
		ParticleBunch(np, generator, beam, filter), GotElastic(false), GotDiffractive(false)
	{
		SetUpProfiling();
	}
	// Proton Bunch Destructor
	//~ProtonBunch(){delete ElasticScatter; delete DiffractiveScatter;};

	virtual bool IsStable() const;
	virtual double GetParticleMass() const;
	virtual double GetParticleMassMeV() const;
	virtual double GetParticleLifetime() const;

	int Scatter(Particle& pi, double x, const Collimator* col);

	int (*ScatterFunctionPointer)(Particle& p, double x, const Collimator* col);

	void set()
	{
		for(int i = 0; i < ntally; tally[i++] = 0)
			;
	}

	void report()
	{
		cout << " Proton Scatter tallies ";
		for(int i = 0; i < ntally; cout << tally[i++] << " ")
			;
		cout << endl;
	}
	/*
	    // set table of t against sigma for calculating b
	    void ConfigureScatter_pp_table(const char*);
	    double get_ft(double t);
	    vector<double> t_sigma_table;
	    double t_sigma_table_step;

	    // set table of xi against sigma for calculating b
	    void ConfigureScatter_xi_table(const char*);
	    double get_fxi(double xi);
	    vector<double> xi_sigma_table;
	    double xi_sigma_table_step;
	 */

	//Scattering physics variables
	double A, Z, E0, X0, rho;
	double lambda_tot;
	double b_pp, b_N;
	double t_low_cut;
	double sigma_pN_total;
	double sigma_pN_elastic;
	double sigma_pn_elastic;
	double sigma_pn_SingleDiffractive;
	double sigma_Rutherford;
	double center_of_mass_squared;
	double I;
	double tmax;
	double C, C0, C1, delta;
	double dEdx;
	double xi0;
	std::string name;
	/**
	 * Select the Scattering physics mode
	 */
	enum scatMode
	{
		SixTrack,
		SixTrackIoniz,
		SixTrackElastic,
		SixTrackSD,
		Merlin

	};

	void EnableScatteringPhysics(scatMode);
	//void EnableSixtrackPhysics(bool);

	int ScatterSixtrack(PSvector& pi, double x, const Collimator* col);
	int ScatterSixtrackAdvancedIonization(PSvector& pi, double x, const Collimator* col);
	int ScatterSixtrackAdvancedElastic(PSvector& pi, double x, const Collimator* col);
	int ScatterSixtrackAdvancedSingleDiffraction(PSvector& pi, double x, const Collimator* col);
	int ScatterMerlin(PSvector& pi, double x, const Collimator* col);

	void ConfigureScatter(const Collimator* col);
	void ConfigureScatterMerlin(const Collimator* col);
	void ConfigureScatterSixtrack(const Collimator* col);
	void ConfigureScatterSixtrackAdvancedIonization(const Collimator* col);
	void ConfigureScatterSixtrackAdvancedElastic(const Collimator* col);
	void ConfigureScatterSixtrackAdvancedSingleDiffraction(const Collimator* col);

	//The pp elastic scattering class
	bool GotElastic;
	ppElasticScatter* ElasticScatter;

	//The pp SingleDiffractive scattering class
	bool GotDiffractive;
	ppDiffractiveScatter* DiffractiveScatter;

private:
	void SetUpProfiling() const;
}; // end ProtonBunch class

} // end namespace ParticleTracking
#endif
