#include "QuadReferenceOutput.h"
#include <fstream>
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/Components.h"
#include <iomanip>

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;


#define WRITE_FOS(w,p,data) (*fosptr)<<scientific<<setw(w)<<setprecision(p)<<(data)

void QuadReferenceOutput::Record(const ComponentFrame* frame, const Bunch* bunch)
{
	if(!frame->IsComponent())
		return;

	const TWRFStructure* cavity = dynamic_cast<const TWRFStructure*>(&(frame->GetComponent()));
	if(cavity) {
		// need to update the reference energy
		refEnergy += (cavity->GetVoltage())/GeV;
		double q = (bunch->GetTotalCharge())*ElectronCharge;
		refEnergy-= eloss*q*Volt;
	}
	else { // must be a quadrupole
		const Quadrupole& quad = static_cast<const Quadrupole&>(frame->GetComponent());
		
		PSvector S;
		bunch->GetCentroid(S);
		double p0 = bunch->GetReferenceMomentum();
		p0*=1+S.dp();

		(*fosptr)<<setw(10)<<left<<quad.GetName()<<right;
		WRITE_FOS(17,8,frame->GetPosition()/1000.0);
		WRITE_FOS(17,8,refEnergy);
		WRITE_FOS(17,8,p0);
		WRITE_FOS(17,8,quad.GetFieldStrength());
		(*fosptr)<<endl;
    }
}

bool QuadReferenceOutput::NewFile(const std::string& fname)
{
	if(fosptr!=0)
		delete fosptr;
	fosptr = new ofstream(fname.c_str());
	return *fosptr;
}
