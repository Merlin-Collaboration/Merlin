#include "TrackingSimulation.h"
#include <fstream>
#include <set>
#include <string>
#include "StringPattern.h"


class TrackingOutput : public SimulationOutput
{
public:
	TrackingOutput(const std::string& fname)
		: SimulationOutput(),fosptr(nullptr)
	{
		NewFile(fname);
	}

	// Close current file and open a new one
	bool NewFile(const std::string& fname);

	// Access the current filestream directly
	ostream& os()
	{
		return *fosptr;
	}

protected:

	virtual void Record(const ComponentFrame* frame, const Bunch* bunch);
	virtual void RecordInitialBunch(const Bunch* bunch)
	{
		Record("INITIAL",bunch);
	}
	virtual void RecordFinalBunch(const Bunch* bunch)
	{
		Record("FINAL",bunch);
	}

private:

	std::ofstream* fosptr;
	void Record(const string&,const Bunch*);
};







