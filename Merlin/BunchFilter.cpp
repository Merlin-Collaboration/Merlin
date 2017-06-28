#include "BunchFilter.h"

namespace ParticleTracking
{


ParticleBunchFilter::~ParticleBunchFilter ()
{
	// Nothing to do
}

bool HorizontalHaloParticleBunchFilter::Apply(const PSvector& v) const
{
	if(v.x() > (orbit+limit) || v.x() < (orbit-limit) )
	{
		if (v.x() != 0.0 && v.xp() != 0.0)
		{
			return true;
		}
	}
	return false;
}

void HorizontalHaloParticleBunchFilter::SetHorizontalLimit(double lim)
{
	cout << "Setting Horizontal limit to: " << lim << endl;
	limit = lim;
}

void HorizontalHaloParticleBunchFilter::SetHorizontalOrbit(double lim)
{
	cout << "Setting Horizontal orbit to: " << lim << endl;
	orbit = lim;
}

bool VerticalHaloParticleBunchFilter::Apply(const PSvector& v) const
{
	if(fabs(v.y()) > limit)
	{
		return true;
	}
	return false;
}

void VerticalHaloParticleBunchFilter::SetVerticalLimit(double lim)
{
	cout << "Setting Vertical limit to: " << lim << endl;
	limit = lim;
}

}
