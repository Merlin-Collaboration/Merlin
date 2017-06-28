#ifndef SymplecticIntegrators_h
#define SymplecticIntegrators_h 1

#include "StdIntegrators.h"
#include "ParticleComponentTracker.h"


namespace ParticleTracking
{
namespace SYMPLECTIC
{

using namespace ParticleTracking;

DECL_SIMPLE_INTG(DriftCI,Drift)
DECL_SIMPLE_INTG(TWRFStructureCI,TWRFStructure)
DECL_SIMPLE_INTG(SWRFStructureCI,SWRFStructure)
DECL_SIMPLE_INTG(RectMultipoleCI,RectMultipole)
DECL_SIMPLE_INTG(MarkerCI,Marker)

// from std intergrators
DECL_SIMPLE_INTG(MonitorCI,Monitor)
DECL_SIMPLE_INTG(SolenoidCI,Solenoid)


class SectorBendCI : public ParticleComponentTracker::Integrator<SectorBend>
{
public:
	void TrackStep(double);
	void TrackEntrance();
	void TrackExit();
};


DECL_INTG_SET(ParticleComponentTracker,StdISet)

}  // end namespace SYMPLECTIC
}
#endif
