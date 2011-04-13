#ifndef ProtonBunch_h
#define ProtonBunch_h 1

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include <iostream>
#include "AcceleratorModel/Aperture.h"
#include "NumericalUtils/PhysicalConstants.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
namespace ParticleTracking
{

class ProtonBunch :public ParticleBunch
{
	static const int ntally=6;
	int tally[ntally];

public:
/*
    //	Constructs a ProtonBunch using the specified momentum,
    //	total charge and the particle array. Note that on exit,
    //	particles is empty.
    ProtonBunch (double P0, double Q, PSvectorArray& particles)
     : ParticleBunch(P0, Q, particles, ProtonMass, ProtonMassMeV, -1.0) {};

    //	Read phase space vectors from specified input stream.
    ProtonBunch (double P0, double Q, std::istream& is)
     : ParticleBunch(P0, Q, is, ProtonMass, ProtonMassMeV, -1.0) {};

    //	Constructs an empty ProtonBunch with the specified
    //	momentum P0 and charge per macro particle Qm (default =
    //	+1).
    ProtonBunch (double P0, double Qm = 1)
     : ParticleBunch(P0, Qm, ProtonMass, ProtonMassMeV, -1.0) {};
*/

    //	Constructs a ProtonBunch using the specified momentum,
    //	total charge and the particle array. Note that on exit,
    //	particles is empty.
    ProtonBunch (double P0, double Q, PSvectorArray& particles)
     : ParticleBunch(P0, Q, particles) {rng();};

    //	Read phase space vectors from specified input stream.
    ProtonBunch (double P0, double Q, std::istream& is)
     : ParticleBunch(P0, Q, is) {rng();};

    //	Constructs an empty ProtonBunch with the specified
    //	momentum P0 and charge per macro particle Qm (default =
    //	+1).
    ProtonBunch (double P0, double Qm = 1)
     : ParticleBunch(P0, Qm) {rng();};

	virtual bool IsStable() const;
	virtual double GetParticleMass() const;
	virtual double GetParticleMassMeV() const;
	virtual double GetParticleLifetime() const;

	int Scatter(PSvector& pi, double x, const Aperture* ap);

	const gsl_rng_type* T;
	gsl_rng* rnd;

	void rng()
	{
		cout << "rng config" << endl;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rnd = gsl_rng_alloc (T);
	}	

	void set()
	{
		for(int i=0;i<ntally;tally[i++]=0);
	}

	void report()
	{
		cout<<" Proton Scatter tallies ";
		for(int i=0; i<ntally; cout << tally[i++] << " ");
		cout<<endl;
	}

    void ConfigureScatter(const Aperture* ap);

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

    //Scattering physics variables
    double A,Z,E0,X0,rho;
    double lambda_tot;
    double b_pp,b_N;
    double t_low_cut;
    double sigma_pN_total;
    double sigma_pN_elastic;
    double sigma_pn_elastic;
    double sigma_pn_SingleDiffractive;
    double sigma_Rutherford;
    double center_of_mass_squared;
    double I;
    double tmax;
    double C,C0,C1,delta;
    double dEdx;
    double xi0;
    


}; // end ProtonBunch class

} // end namespace ParticleTracking
#endif
