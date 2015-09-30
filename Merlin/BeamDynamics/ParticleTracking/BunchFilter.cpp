#include "BeamDynamics/ParticleTracking/BunchFilter.h"

namespace ParticleTracking {


ParticleBunchFilter::~ParticleBunchFilter ()
{
    // Nothing to do
}

bool HorizontalHaloParticleBunchFilter::Apply(const PSvector& v) const
{
	if(v.x() > (orbit+limit) || v.x() < (orbit-limit) )
	{
		if (v.x() != 0.0 && v.xp() != 0.0){
			//cout << v.x() << "\t" << limit << endl;
			return true;
		}
	}
	//cout << v.x() << "\t" << limit << endl;
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
//		cout << "Particle passes filter" << endl;
//		cout << v.y() << "\t" << limit << endl << endl;
		return true;
	}
//	cout << "Particle fails filter" << endl;
//	cout << v.y() << "\t" << limit << endl << endl;
	return false;
}

void VerticalHaloParticleBunchFilter::SetVerticalLimit(double lim)
{
	cout << "Setting Vertical limit to: " << lim << endl;
	limit = lim;
}

}
