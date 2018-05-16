#include <cmath>
#include "ElectronBunch.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "MerlinException.h"
#include "Aperture.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

double ElectronBunch::GetParticleMass() const
{
	return ElectronMass;
}

double ElectronBunch::GetParticleMassMeV() const
{
	return ElectronMassMeV;
}

double ElectronBunch::GetParticleLifetime() const
{
	return 0;
}

bool ElectronBunch::IsStable() const
{
	return true;
}
/*
int ElectronBunch::Scatter(PSvector& p,double x,double E0,const Aperture* ap)
{
	return 0;
} //End of ScatterElectron
*/
