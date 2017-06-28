#ifndef FlukaCollimationOutput_h
#define FlukaCollimationOutput_h 1

#include <string>
#include <vector>

#include "CollimationOutput.h"
#include "AcceleratorComponent.h"
#include "ParticleBunch.h"
#include "PSTypes.h"

namespace ParticleTracking
{

class FlukaCollimationOutput : public CollimationOutput
{

public:

	FlukaCollimationOutput(OutputType otype = tencm);
	~FlukaCollimationOutput();

	virtual void Finalise();
	virtual void Output(std::ostream* os);
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

protected:

private:

};

}

#endif

