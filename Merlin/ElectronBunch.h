#ifndef ElectronBunch_h
#define ElectronBunch_h 1

#include "ParticleBunch.h"
#include <iostream>
#include "Aperture.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
namespace ParticleTracking
{

class ElectronBunch :public ParticleBunch
{
	static const int ntally=6;
	int tally[ntally];

public:

	/**
	*	Constructs an ElectronBunch using the specified momentum,
	*	total charge and the particle array. Note that on exit,
	*	particles is empty.
	*/
	ElectronBunch (double P0, double Q, PSvectorArray& particles)
		: ParticleBunch(P0, Q, particles) {};

	/**
	*	Read phase space vectors from specified input stream.
	*/
	ElectronBunch (double P0, double Q, std::istream& is)
		: ParticleBunch(P0, Q, is) {};

	/**
	*	Constructs an empty ElectronBunch with the specified
	*	momentum P0 and charge per macro particle Qm (default =
	*	+1).
	*/
	ElectronBunch (double P0, double Qm = 1)
		: ParticleBunch(P0, Qm) {};

	virtual bool IsStable() const;
	virtual double GetParticleMass() const;
	virtual double GetParticleMassMeV() const;
	virtual double GetParticleLifetime() const;

	//int Scatter(PSvector& pi, double x, double E0, const Aperture* ap);
	//int Scatter(PSvector&, double, double, const Aperture*);

	void set()
	{
		for(int i=0; i<ntally; tally[i++]=0) {}
	}

	void report()
	{
		cout << "Electron Scatter tallies ";
		for(int i=0; i<ntally; cout << tally[i++] << " ") {}
		cout<<endl;
	}
}; // end ElectronBunch class
} // end namespace ParticleTracking
#endif
