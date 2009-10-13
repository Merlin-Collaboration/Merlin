/////////////////////////////////////////////////////////////////////////
// File ModelConstruction.h implementation
// Routine to construct ILC Main Linac model
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include <vector>
#include <sstream>

#include "merlin_config.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "XTFFInterface_1.h"
#include "BeamModel/BeamData.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "TeslaWakePotential.h"
#include "AcceleratorModel/Frames/PatchFrame.h"
#include "AcceleratorModel/ControlElements/Klystron.h"
#include "AcceleratorModel/Frames/TComponentFrame.h"
#include "ILCDFS_IO.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname, bool addcurv)
{
	// Charge per bunch needed to estimate beam loading
	const double Qb = 2.0e+10;

	dfs_trace(dfs_trace::level_1)<<"Constructing model from "<<fname;
	if(addcurv)
		dfs_trace(dfs_trace::level_1)<<" with Earth curvature";
	dfs_trace(dfs_trace::level_1)<<endl;


	ofstream logs("construction.log");
	XTFFInterface_1 mc(fname,Qb,&logs);
	mc.ConstructGirders(true);

	pair<AcceleratorModel*,BeamData*> mb = mc.Parse();
	BeamData* beam0 = mb.second;
	AcceleratorModel* model = mb.first;

	// The following quantities are not
	// int the TAPE file
	double gamma = beam0->p0/MeV/ElectronMassMeV;
	beam0->emit_x = 8.0e-06/gamma;
	beam0->emit_y = 0.02e-06/gamma;
	beam0->charge = Qb;
	beam0->sig_dp = 0.0107; // from PT's bunch compressor
	beam0->sig_z  = 300.0e-06;

	// Now add wakefields to cavities

	vector<TComponentFrame<TWRFStructure>*> cavities;
	model->ExtractTypedComponents(cavities,"*"); // linac cavities only

	TeslaWakePotentials* wake = new TeslaWakePotentials;	
	CircularAperture* iris = new CircularAperture(0.035); // TESLA iris aperture

	size_t nCavsPerKlystron = 24;
	size_t nKlys=0;
	vector<RFStructure*> kcavs;

	for(size_t n=0; n<cavities.size(); n++) {

		TWRFStructure& cav = cavities[n]->GetComponent();
		cav.SetWakePotentials(wake);
		cav.SetAperture(iris);

		kcavs.push_back(&cav);
		if(kcavs.size()==nCavsPerKlystron) {
			nKlys++;
			ostringstream ss;
			ss<<"KLYS_"<<nKlys;
			model->AddModelElement(new Klystron(ss.str(),kcavs));
			kcavs.clear();
		}
	}

	// Add last klystron
	if(!kcavs.empty()) {
		nKlys++;
		ostringstream ss;
		ss<<"KLYS_"<<nKlys;
		model->AddModelElement(new Klystron(ss.str(),kcavs));
		kcavs.clear();
	}

	if(addcurv) {

		// Now adjust the  machine to approximately follow the earth's curvature
		// We do this by extracting the VPIV geometry patch frames.

		AcceleratorModel::Beamline let = model->GetBeamline();
		PatchFrame* vpatch0 = 0;
		double earthRho = 6.0e+06; // 6000 km
		double totalAngle = 0.0;
		double totalDist = 0.0;
		size_t npatches = 0;

		for(AcceleratorModel::BeamlineIterator ci = let.begin(); ci!=let.end(); ci++) {

			PatchFrame* vpatch = dynamic_cast<PatchFrame*>(*ci);

			if(vpatch) {
				npatches++;

				if(vpatch0==0)
					vpatch0 = vpatch; // first vpatch
				else {

					double dist = (*vpatch).GetPosition() - (*vpatch0).GetPosition();
					double angle = dist/earthRho;
					GeometryPatch* gp = new GeometryPatch;
					gp->RotateX(angle);
					vpatch->SetGeometryPatch(gp); // rotates down
					totalAngle+=angle;
					totalDist +=dist;
					vpatch0 = vpatch;
				}
			}
		}
		//		cout<<"number of vertical patches: "<<npatches<<endl<<endl;
		//		cout<<"total angle:  "<<totalAngle/degree<<endl;
		//		cout<<"total length: "<<totalDist<<endl;
		//		cout<<"radius:       "<<totalDist/totalAngle/1000<<" km\n"<<endl;

		// Approximate match beam conditions
		beam0->Dy = 0.0007;
		beam0->Dyp = -0.00001;
	}
	return mb;
}

