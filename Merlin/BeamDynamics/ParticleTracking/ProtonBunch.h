#ifndef ProtonBunch_h
#define ProtonBunch_h 1

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include <iostream>
#include "AcceleratorModel/Aperture.h"
#include "NumericalUtils/PhysicalConstants.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
namespace ParticleTracking
{

class ProtonBunch :public ParticleBunch
{
	static const int ntally=6;
	int tally[ntally];

    //	Constructs a ProtonBunch using the specified momentum,
    //	total charge and the particle array. Note that on exit,
    //	particles is empty.
    ProtonBunch (double P0, double Q, PSvectorArray& particles, const double ParticleMass = ProtonMass, const double ParticleMassMev = ProtonMassMeV, const double ParticleLifetime = -1)
     : ParticleBunch(P0, Q, particles) {};

    //	Read phase space vectors from specified input stream.
    ProtonBunch (double P0, double Q, std::istream& is, const double ParticleMass = ProtonMass, const double ParticleMassMev = ProtonMassMeV, const double ParticleLifetime = -1)
     : ParticleBunch(P0, Q, is) {};

    //	Constructs an empty ProtonBunch with the specified
    //	momentum P0 and charge per macro particle Qm (default =
    //	+1).
    ProtonBunch (double P0, double Qm = 1, const double ParticleMass = ProtonMass, const double ParticleMassMev = ProtonMassMeV, const double ParticleLifetime = -1)
     : ParticleBunch(P0, Qm) {};

	virtual bool IsStable() const;

public:

	int Scatter(PSvector& pi, double x, double E0, const Aperture* ap);

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
}; // end ProtonBunch class
}; // end namespace ParticleTracking
#endif
