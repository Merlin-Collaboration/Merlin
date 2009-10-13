#include "NumericalUtils/PhysicalUnits.h"
#include "MADInterface/MADInterface.h"
#include "AcceleratorModel/StdComponent/StandardMultipoles.h"
#include "AcceleratorModel/Supports/MagnetMover.h"

#include "RingDynamics/LatticeFunctions.h"

#define BEAMENERGY 5.0*GeV

typedef vector<MagnetMover*> MagnetMoverList;
typedef vector<Quadrupole*> QuadList;

using namespace PhysicalUnits;

int main()
{
	// Construct the AcceleratorModel
	// from a lattice file produced by MAD
	MADInterface madi("../lattices/MERLINFodo.lattice.txt", BEAMENERGY);

	ofstream madlog("mad.log");
	madi.SetLogFile(madlog);
	madi.SetLoggingOn();

	AcceleratorModel* theModel = madi.ConstructModel();


	// Extract a list of magnet movers from the AcceleratorModel
	// and translate the 20th mover 20 microns vertically.
	MagnetMoverList magnetMovers;
	theModel->ExtractTypedElements(magnetMovers);
	magnetMovers[20]->SetY(20.0e-6);


	// Extract a list of the quadrupoles,
	// and put a 5% gradient error on the 20th quadrupole.
	QuadList quads;
	theModel->ExtractTypedElements(quads);
	MultipoleField& field = quads[20]->GetField();
	Complex b1 = field.GetComponent(1);
	field.SetComponent(1, b1.real()*1.05, b1.imag()*1.05);


	// Find the lattice functions.
	// Note the default is to find the closed orbit first,
	// and use the transfer matrix around the closed orbit.
	// Default functions are the closed orbit co-ordinates
	// and the coupled-lattice equivalents of the 
	// Twiss alpha and beta functions.
	LatticeFunctionTable latticeFunctions(theModel,BEAMENERGY);

	// We add functions to give us the dispersion.
	// Dx  = beta(1,6,3)/beta(6,6,3)
	// Dpx = beta(2,6,3)/beta(6,6,3)
	// etc.
	latticeFunctions.AddFunction(1,6,3);
	latticeFunctions.AddFunction(2,6,3);
	latticeFunctions.AddFunction(3,6,3);
	latticeFunctions.AddFunction(4,6,3);
	latticeFunctions.AddFunction(6,6,3);

	// Calculate the lattice functions.
	latticeFunctions.Calculate();

	// Generate a file with the results.
	// First column is the s position in the beamline.
	ofstream latticeFunctionLog("LatticeFunctions.dat");
	latticeFunctions.PrintTable(latticeFunctionLog);
	
	delete theModel;

	cout<<"Finished!"<<endl;
	return 0;
}
