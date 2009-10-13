#include "BeamDynamics/TrackingSimulation.h"
#include <fstream>
#include <string>


class QuadReferenceOutput : public SimulationOutput {
public:
	QuadReferenceOutput(const std::string& fname, double E0, double el) 
		: SimulationOutput(),fosptr(0),refEnergy(E0),eloss(el) {
			NewFile(fname);
			AddIdentifier("Quadrupole.*");
			AddIdentifier("TWRFStructure.*");
			output_all=output_initial=output_final=false;
	}

	// Close current file and open a new one
	bool NewFile(const std::string& fname);

	// Access the current filestream directly
	ostream& os() { return *fosptr; }

protected:

	virtual void Record(const ComponentFrame* frame, const Bunch* bunch);
	virtual void RecordInitialBunch(const Bunch* bunch) {}
	virtual void RecordFinalBunch(const Bunch* bunch){}

private:

	std::ofstream* fosptr;
	double refEnergy;
	double eloss;
};







