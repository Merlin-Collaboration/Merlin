#include <cmath>
#include "MuonBunch.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "MerlinException.h"
#include "Aperture.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

double MuonBunch::GetParticleMass() const
{
	return MuonMass;
}

double MuonBunch::GetParticleMassMeV() const
{
	return MuonMassMeV;
}

double MuonBunch::GetParticleLifetime() const
{
	return MuonLifetime;
}

bool MuonBunch::IsStable() const
{
	return false;
}
/*
int MuonBunch::Scatter(PSvector& p,double x,double E0,const Aperture* ap)
{
	return 0;
} //End of ScatterMuon
*/
