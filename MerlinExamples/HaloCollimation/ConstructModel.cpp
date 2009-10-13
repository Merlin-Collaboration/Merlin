// ConstructModel.cpp
// --------------------------------------------------------------------
// Implementation of global function to construct
// an accelerator model from an OPTICS MAD output file.
//

#include "merlin_config.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "MADInterface/MADInterface.h"
#include "BeamModel/BeamData.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/Miscellaneous/CorrectorWinding.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname)
{
    // Here we
    // 1. construct the initial (matched) beam data
    // 2. construct the AcceleratorModel from the supplied
    //    file.

    double P0 = 250.0*GeV;
	double gamma = P0/(ElectronMassMeV*MeV);

	BeamData *beam = new BeamData;
	beam->beta_x = 172.070;
	beam->sig_z  = 0.0003;
	beam->beta_y = 57.48;
	beam->sig_dp = 1.5e-03;
	beam->emit_x = 10.0e-06/gamma;
	beam->emit_y = 0.03e-06/gamma;
	beam->charge = 2.0e+10;
	beam->p0=P0;

	MADInterface mad(fname,P0);	

	ofstream madlog("mad.log");
	mad.SetLogFile(madlog);
	mad.SetLoggingOn();

	mad.ConstructApertures(true);
	mad.ConstructFlatLattice(true);

    // MADInterface knows nothing about the following types,
    // so we treat them as drifts.
	mad.TreatTypeAsDrift("ABSORBER");
	mad.TreatTypeAsDrift("SAMPLER");
	mad.TreatTypeAsDrift("SPOILER");

	AcceleratorModel* model = mad.ConstructModel();
	return make_pair(model,beam);
}
