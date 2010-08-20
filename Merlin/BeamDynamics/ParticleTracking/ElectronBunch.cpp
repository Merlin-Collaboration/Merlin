#include <cmath>
#include "ElectronBunch.h"
//#include "Collimators/CoulombScatter.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "Exception/MerlinException.h"
#include "AcceleratorModel/Aperture.h"
#include "AcceleratorModel/Apertures/TiltedAperture.hpp"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

bool ElectronBunch::IsStable() const
{
	return true;
}

int ElectronBunch::Scatter(PSvector& p,double x,double E0,const Aperture* ap)
{
	return 0;
} //End of ScatterElectron
