/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
using namespace std;

/*
    Tests materials class
      and also use of the CERN ROOT histogram package

      In order to do this, you need to set the environment variable
      ROOTSYS to the directory where root has been installed
      e.g.  /home/software/root

      This can be cleanly done by sourcing the file bin/thisroot.sh or .csh
      in that directory
 */

#include "MaterialProperties.h"
#include "MaterialData.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

//
#include "BeamData.h"
#include "PSmoments.h"
#include "ParticleBunch.h"
#include "ParticleDistributionGenerator.h"
#include "ParticleTracker.h"
#include "RandomNG.h"
#include "AcceleratorModel.h"
// #include "SimpleApertures.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "AcceleratorModelConstructor.h"
#include "Drift.h"
#include "CollimateProtonProcess.h"
#include "Collimator.h"
#include "CollimatorAperture.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <string>
#include <map>
using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

TH1D *histt1;
TH1D *histt2;
TH1D *histt3;
TH1D *histt4;
TH1D *histt5;

int main(int argc, char* argv[])
{
	cout << "Here we go ! \n";
   const double Z_pp = 35.4548e-3;
                const double B_pp = 0.30810e-3;
                const double Y1_pp = 42.5323e-3;
                const double Y2_pp = 33.3433e-3;
                const double eta1 = 0.45817;
                const double eta2 = 0.545;
                const double s0 = 28.998225 * PhysicalUnits::GeV * PhysicalUnits::GeV;
                const double s1 = 1 * PhysicalUnits::GeV * PhysicalUnits::GeV;
                const double s = 7000*7000-2*7000*.938;
             cout<<" total cross section calcultion "<<
                Z_pp + B_pp * pow(log(s / s0), 2.0) + Y1_pp * pow(s1 / s, eta1) 
- Y2_pp * pow(s1 / s, eta2)<<endl;;
	int seed = 0;
	try{
		if(argc >= 2)
			seed = atoi(argv[1]);
		cout << " seed set to " << seed << endl;
		RandomNG::init(seed);         // initialise Random number generator
		TFile hfile("job1.root", "RECREATE", "merlin output");
		int npart = 100000;
		if(argc >= 3)
			npart = atoi(argv[2]);
		cout << npart << " particles \n";

		BeamData mybeam;
		mybeam.charge = 1.31e11;
		mybeam.beta_x = 0.5495121695 * meter;
		mybeam.beta_y = 0.5498820579 * meter;

		mybeam.emit_x = 33.640 * 5.026457122e-10 * meter;
		mybeam.emit_y = 33.64 * 5.026457122e-10 * meter;
		mybeam.emit_x = 0 * meter;
		mybeam.emit_y = 0 * meter;
		mybeam.sig_z = 75.5 * millimeter;
		// mybeam.sig_dp = 0.000113 ;
		mybeam.sig_dp = 0.0;
		mybeam.p0 = 7000 * GeV;
		//mybeam.alpha_x = -0.001*meter;
		//mybeam.alpha_y = -0.001*meter;

		mybeam.alpha_x = -0.0001721885021 * meter;
		mybeam.alpha_y = -0.0004654580947 * meter;

		int i;
		mybeam.yp0 = 0;
		mybeam.xp0 = 0;
		mybeam.x0 = 0;
		MaterialProperties xx(1., 2., 3., 4., 5., 6., 7., 8.);
		MaterialProperties* yy =
			new  MaterialProperties(1, 2, 3, 4, 5, 6, 7, 8);
		cout << " Simple Material properties " << (*yy) << endl;
		(*yy->extra)["stuff"] = 99;
		(*yy->extra)["more"] = 98;
		yy->SetExtra("new cheese sandwich", double(88.0), double(99), 55.6);
		cout << "with extras " << (*yy) << endl;
		cout << " test GetExtra and get " << yy->GetExtra("more") << endl;
		StandardMaterialData*  matter = new StandardMaterialData;

		cout << "Standard  materialdata" << endl;
		matter->PrintTable();
///		matter->UseSixTrackValues();
///		cout << "SixTrack Modified  materialdata" << endl;
		matter->PrintTable();
		MaterialProperties test1 = *(matter->property[string("Cu")]);
		cout << " for example copper is " << test1 << endl;
		MaterialProperties test2(*(matter->property["Al"]));
		cout << " and aluminium is " << test2 << endl;
		MaterialProperties test3;
		test3 = test2;
		cout << "test  copy of copper and get  " << test3 << endl;
		(*(matter->property["Cu"])->extra)["slope"] = 4.;
		cout << " now add property to copper\n  ";
		matter->PrintTable();
		matter->MakeMixture("SiC", "Si C", 1., 1., 3.2 * gram / cc);
		matter->MakeMixtureByWeight("CuDiamond", "Cu C", .63, .37, 3.2 * gram / cc);
		matter->MakeMixture("test", "Al Be Cu W", 1., 2., 3., 4., 99., 88.);
		matter->MakeMixtureByWeight("test2", "Al Be Cu W", 1., 2., 3., 4., 99., 88.);
		matter->PrintTable();
		Aperture* ap = new CircularAperture(.2);
		ParticleDistributionGenerator* pg = new NormalParticleDistributionGenerator();

		double lim1[3] = {0.00001, 0.0001, 0.0001};
		double lim2[3] = {0.00001, 0.0001, 0.0001};
		double lim3[3] = {0.00001, 0.0005, 0.0001};

		string trymaterial[3] = {"Cu", "C", "CuDiamond"};
		string types[2] = {"Full", "Half"};
		double thickness[] = {1., 10., 30.};
// loop over type
		for(int itype = 0; itype < 2; itype++)
		{
			cout << "Aperture " << types[itype] << endl;
// loop over thickness
			for(int ithick = 0; ithick < 3; ithick++)
			{
// loop over material
				for(int imat = 0; imat < 3; imat++)
				{
					double offset = (itype == 0) ? 0 : (1.0 + 1.E-6) * meter;
					mybeam.y0 = offset;
					ParticleBunch* myBunch = new ParticleBunch(npart, NormalParticleDistributionGenerator(), mybeam);
					cout << " material " << trymaterial[imat] << endl;
					string histname = types[itype] + "." + to_string(int(thickness[ithick])) + "." + trymaterial[imat];
					cout << " hist name " << histname << endl;
					cout << " y offset " << offset << endl;
					string tempstring;
					tempstring = histname + ".xbefore";
					TH1D *PShist1 = new TH1D(tempstring.c_str(), "x before", 100, -lim1[ithick], lim1[ithick]);
					tempstring = histname + ".xafter";
					TH1D *PShist2 = new TH1D(tempstring.c_str(), "x after", 100, -lim1[ithick], lim2[ithick]);
					tempstring = histname + ".xpbefore";
					TH1D *PShist3 = new TH1D(tempstring.c_str(), "x prime before", 100, -lim2[ithick], lim2[ithick]);
					tempstring = histname + ".xpafter";
					TH1D *PShist4 = new TH1D(tempstring.c_str(), "x prime after", 100, -lim2[ithick], lim2[ithick]);
					tempstring = histname + ".dpafter";
					TH1D *PShist5 = new TH1D(tempstring.c_str(), "delta p after", 100, 0, lim3[ithick]);
					tempstring = histname + ".yafter";
					TH1D *PShist6 = new TH1D(tempstring.c_str(), "y after", 100, -lim1[ithick], lim1[ithick]);
					tempstring = histname + ".ypafter";
					TH1D *PShist7 = new TH1D(tempstring.c_str(), "y prime after", 100, -lim2[ithick], lim2[ithick]);

					PSvectorArray particlearray1 = myBunch->GetParticles();

					// int nafter = particlearray1.size();

					for(int i = 0; i <= npart; i++)
					{
						PShist3->Fill(particlearray1[i].xp());
						PShist1->Fill(particlearray1[i].x());
					}

					AcceleratorModelConstructor* myaccmodelctor = new AcceleratorModelConstructor()
					;
					myaccmodelctor->NewModel();
					myaccmodelctor->AppendComponent(new Drift("DRIFT1", 1.0 * meter));
					Collimator* TestCol = new Collimator("TheCollimator", thickness[ithick] * 0.01);
					TestCol->SetMaterialProperties(matter->property[trymaterial[imat]]);
					CollimatorAperture* app;
					if(itype == 0)
					{
						app = new CollimatorAperture(0, 0, 0, .01, 0, 0);
					}
					else
					{
						app = new CollimatorAperture(2, 2, 0, .01, 0, 0);
					}
					//if(itype==0){
					//app->SetExitWidth(0);
					//app->SetExitHeight(0);
					//app->SetEntranceWidth(0);
					//app->SetEntranceHeight(0);
					//} else {
					//app->SetExitWidth(2);
					//app->SetExitHeight(2);
					// }
					myaccmodelctor->AppendComponent(TestCol, 0.01);
					myaccmodelctor->AppendComponent(new Drift("DRIFT1", 1.0 * meter));
					TestCol->SetAperture(app);
					AcceleratorModel* mymodel = myaccmodelctor->GetModel();
					ParticleTracker mytracker(mymodel->GetBeamline(), myBunch);
					ParticleTracker* tracker = new ParticleTracker(mymodel->GetBeamline(), myBunch
						);

					CollimateProtonProcess* myCollimateProcess = new CollimateProtonProcess(0, 7);
					ScatteringModel* myScatter = new ScatteringModel();
					myCollimateProcess->SetScatteringModel(myScatter);
					myCollimateProcess->ScatterAtCollimator(true); // Needs resurrection
					tracker->AddProcess(myCollimateProcess);
					ParticleBunch* bunch2;
					cout << " start tracking" << endl;
					bunch2 = tracker->Track(myBunch);
					cout << " done tracking" << endl;
					PSvectorArray myparticles2 = bunch2->GetParticles();
					int nafter = myparticles2.size();
					cout << " number surviving " << nafter << endl;
					for(i = 0; i <= nafter; i++)
					{
						PShist2->Fill(myparticles2[i].x());
						PShist4->Fill(myparticles2[i].xp());
						PShist5->Fill(-myparticles2[i].dp());
						PShist6->Fill(myparticles2[i].y() - offset);
						PShist7->Fill(myparticles2[i].yp());
					}
				}
			}
		}
		cout << "writing\n";
		hfile.Write();
		cout << "written\n";
		delete matter;
		cout << "deleted\n";
	}     // end of try
	catch(MerlinException&  s)
	{
		cout << " Merlin Exception: " << s.Msg() << endl;
	}

	return 0;

}
