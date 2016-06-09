#ifndef LossMapCollimationOutput_h
#define LossMapCollimationOutput_h 1

#include <string>
#include <vector>

#include "Collimators/Output/CollimationOutput.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamModel/PSTypes.h"
namespace ParticleTracking
{

class LossMapCollimationOutput : public CollimationOutput
{

public:

	LossMapCollimationOutput(OutputType otype = tencm);
	~LossMapCollimationOutput();

	virtual void Finalise();
	virtual void Output(std::ostream* os);
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

	void SetWarmRegion(std::pair<double, double>);
	void ClearWarmRegions();
	std::vector<std::pair<double,double> > GetWarmRegions() const;

protected:
	std::vector<std::pair<double,double> > WarmRegions;

private:

};

}
#endif

