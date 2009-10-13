#include "merlin_config.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "MADInterface/XTFFInterface.h"
#include "BeamModel/BeamData.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "TeslaWakePotential.h"
#include "AcceleratorModel/Frames/PatchFrame.h"
#include <vector>

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname)
{
	double qt = 2.0e+10; // single-bunch charge

	ofstream logs("construction.log");
	XTFFInterface mc(fname,qt,&logs);
	mc.ConstructGirders(true);

//	mc.TreatTypeAsDrift("LCAV");
	
	pair<AcceleratorModel*,BeamData*> mb = mc.Parse();
	BeamData* beam0 = mb.second;
	AcceleratorModel* model = mb.first;
	
	// The following quantities are not
	// int the TAPE file
	
	double gamma = beam0->p0/MeV/ElectronMassMeV;
	beam0->emit_x = 8.0e-06/gamma;
	beam0->emit_y = 0.02e-06/gamma;
	beam0->charge = qt==0 ? 1.0 : qt;
	beam0->sig_dp = 0.028;
	beam0->sig_z  = 300.0e-06;
	
	// Now add wakefields to cavities
	
	vector<TWRFStructure*> cavities;
	model->ExtractTypedElements(cavities,"CAV*"); // linac cavities only
	TeslaWakePotentials* wake = new TeslaWakePotentials;
	
	CircularAperture* iris = new CircularAperture(0.035); // TESLA iris aperture
	
	for(vector<TWRFStructure*>::iterator c = cavities.begin(); c!=cavities.end(); c++) {
		(*c)->SetWakePotentials(wake);
		(*c)->SetAperture(iris);
	}
	
	return mb;
}

