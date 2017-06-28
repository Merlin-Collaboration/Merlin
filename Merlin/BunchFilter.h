#ifndef _BUNCHFILTER_H_
#define _BUNCHFILTER_H_

//#include "ParticleBunchConstructor.h"
#include "PSvector.h"

namespace ParticleTracking
{

class ParticleBunchFilter
{
public:

	virtual ~ParticleBunchFilter ();

	//	Used by a ParticleBunchConstructor object to select
	//	vectors for inclusion in a ParticleBunch.
	virtual bool Apply (const PSvector& v) const = 0;
};

class HorizontalHaloParticleBunchFilter : public ParticleBunchFilter
{
public:

	//	Used by a ParticleBunchConstructor object to select
	//	vectors for inclusion in a ParticleBunch.
	bool Apply (const PSvector& v) const;

	void SetHorizontalLimit(double);
	void SetHorizontalOrbit(double);

private:
	double limit;
	double orbit;
};


class VerticalHaloParticleBunchFilter : public ParticleBunchFilter
{
public:

	//	Used by a ParticleBunchConstructor object to select
	//	vectors for inclusion in a ParticleBunch.
	bool Apply (const PSvector& v) const;

	void SetVerticalLimit(double);

private:
	double limit;
};


}	//End particle tracking namespace

#endif
