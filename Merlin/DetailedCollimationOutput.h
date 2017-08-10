#ifndef DetailedCollimationOutput_h
#define DetailedCollimationOutput_h 1

#include <string>
#include <vector>

#include "CollimationOutput.h"
#include "AcceleratorComponent.h"
#include "PSTypes.h"
#include "StringPattern.h"

namespace ParticleTracking
{

class DetailedCollimationOutput : public CollimationOutput
{
public:
	DetailedCollimationOutput(OutputType otype = tencm);
	~DetailedCollimationOutput();

	virtual void Finalise();
	virtual void Output(std::ostream* os);
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

	virtual void AddIdentifier(const std::string e);

private:
	std::vector<StringPattern> ids;

};

}

#endif

