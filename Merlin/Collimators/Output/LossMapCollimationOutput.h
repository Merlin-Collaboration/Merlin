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

	/**
	* Finalise will call any sorting algorithms and perform formatting for final output
	*/
	virtual void Finalise();

	/**
	* Outputs the loss map data to a specified output stream
	* @param[out] os The stream to output to.
	*/
	virtual void Output(std::ostream* os);

	/**
	* Called from CollimateProtonProcess::DeathReport to add a particle to the CollimationOutput
	*/
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

	/**
	* Sets a warm area of the machine.
	* @param[in] wr A std::pair that contains the start and end location of a warm region. First contains the start location, and second the end.
	*/
	void SetWarmRegion(std::pair<double, double> wr);

	/**
	* Clears out any previously defined warm regions.
	*/
	void ClearWarmRegions();

	/**
	* Gets the vector containing the currently set warm regions of the machine.
	* @return The std::vector of std::pair<double,double> with each entry containing the start and end locations of warm regions.
	*/
	std::vector<std::pair<double,double> > GetWarmRegions() const;

protected:

	//A vector of std::pair containing the start and end of warm regions of the machine. Can be empty. First contains the start location, and second the end.
	std::vector<std::pair<double,double> > WarmRegions;

private:

};

}
#endif

