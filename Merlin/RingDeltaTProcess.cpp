#include <cmath>
#include <algorithm>
#include "utils.h"

#include "RingDeltaTProcess.h"
#include "SectorBend.h"

namespace ParticleTracking
{

struct ApplyDeltaT
{
private:
	double dt;

public:
	ApplyDeltaT(double _dt) : dt(_dt) {};

	void operator()(PSvector& v)
	{
		v.ct() += dt;
	};
};

RingDeltaTProcess::RingDeltaTProcess (int prio)
	: ParticleBunchProcess("RING DELTA T",prio) {}

void RingDeltaTProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	active = (dynamic_cast<SectorBend*>(&component)) ? true : false;
	dL     = component.GetLength();
	intS   = 0;
}

void RingDeltaTProcess::DoProcess (double ds)
{
	intS += ds;
	for_each(currentBunch->begin(),currentBunch->end(),ApplyDeltaT(scale*ds));
	active = intS!=dL;
}

double RingDeltaTProcess::GetMaxAllowedStepSize () const
{
	return dL-intS;
}

void RingDeltaTProcess::SetBendScale (double bendscale)
{
	scale = bendscale;
}

} // end namespace ParticleTracking
