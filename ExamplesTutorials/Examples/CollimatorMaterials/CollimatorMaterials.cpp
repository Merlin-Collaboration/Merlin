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
#include "CollimateParticleProcess.h"

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
		// mybeam.y0 = offset;

		int i;
		float offset = (2.10 + 1.E-6) * meter;
		mybeam.yp0 = 0;
		mybeam.xp0 = 0;
		mybeam.x0 = 0;
		mybeam.y0 = offset;
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
		matter->UseSixTrackValues();
		cout << "SixTrack Modified  materialdata" << endl;
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
		matter->MakeMixture("SiC", "Si C", 1, 1, 3.2, 4.3);
		cout << " done 1" << endl;
		matter->MakeMixture("test", "Al Be Cu W", 1, 2, 3, 4, 99., 88.);
		matter->MakeMixtureByWeight("test2", "Al Be Cu W", 1, 2, 3, 4, 99., 88.);
		matter->PrintTable();
		Aperture* ap = new CircularAperture(.2);
		ParticleDistributionGenerator* pg = new NormalParticleDistributionGenerator();

		ParticleBunch* myBunch = new ParticleBunch(npart, NormalParticleDistributionGenerator(), mybeam);

		double xlim, ylim, xplim, yplim, zlim;
		xlim  = 0.0001;
		xplim = 0.0001;
		ylim  = sqrt(mybeam.emit_y * mybeam.beta_y) * 4.0;
		yplim = sqrt(mybeam.emit_y / mybeam.beta_y) * 4.0;
		ylim = 0.00001;
		yplim = 0.00003;
		zlim = mybeam.sig_z * 4.0;

		histt1 = new TH1D("t1", "Nuclear elastic t", 100, 0, 0.2);
		histt2 = new TH1D("t2", "Nucleon elastic t", 100, 0, 0.2);
		histt3 = new TH1D("t3", "SD Mass squared", 100, 0, 100.0);
		histt4 = new TH1D("t4", "Diffractive  t", 100, 0, 0.2);
		histt5 = new TH1D("m1", "Diffractive Mass squared", 100, 0, 50);
		TH1D *PShist1 = new TH1D("xbefore", "x before", 100, -xlim, xlim);
		TH1D *PShist2 = new TH1D("xafter", "x after", 100, -xlim, xlim);
		TH1D *PShist3 = new TH1D("xpbefore", "x prime before", 100, -xplim, xplim);
		TH1D *PShist4 = new TH1D("xpafter", "x prime after", 100, -xplim, xplim);
		TH1D *PShist5 = new TH1D("dpafter", "delta p after", 100, 0, 0.001);
		TH1D *PShist6 = new TH1D("yafter", "y after", 100, -xlim, xlim);
		TH1D *PShist7 = new TH1D("ypafter", "y prime after", 100, -xplim, xplim);
		//   TH2D *yPShist3 = new TH2D("test1","test2", 100, -zlim,zlim , 100, -yplim, yplim);

		PSvectorArray particlearray1 = myBunch->GetParticles();

		npart = particlearray1.size();

		for(int i = 0; i <= npart; i++)
		{
			PShist3->Fill(particlearray1[i].xp());
			PShist1->Fill(particlearray1[i].x());
		}

		AcceleratorModelConstructor* myaccmodelctor = new AcceleratorModelConstructor()
		;
		myaccmodelctor->NewModel();
		myaccmodelctor->AppendComponent(new Drift("DRIFT1", 1.0 * meter));
		AcceleratorModel* mymodel = myaccmodelctor->GetModel();
		ParticleTracker mytracker(mymodel->GetBeamline(), myBunch);
		ParticleTracker* tracker = new ParticleTracker(mymodel->GetBeamline(), myBunch
			);

		CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(0, 7);
		myCollimateProcess->ScatterAtCollimator(true);          // Needs resurrection
		tracker->AddProcess(myCollimateProcess);
		ParticleBunch* bunch2;
		cout << " start tracking" << endl;
		bunch2 = tracker->Track(myBunch);
		cout << " done tracking" << endl;
		cout << " get particles \n";
		PSvectorArray myparticles2 = bunch2->GetParticles();
		npart = myparticles2.size();
		cout << " size " << npart << endl;
		for(i = 0; i <= npart; i++)
		{
			PShist2->Fill(myparticles2[i].x());
			PShist4->Fill(myparticles2[i].xp());
			PShist5->Fill(-myparticles2[i].dp());
			PShist6->Fill(myparticles2[i].y() - offset);
			PShist7->Fill(myparticles2[i].yp());
////                        PShist3->Fill((myparticles2[i]).ct(),myparticles2[i].yp());

		}
		PShist2->Draw();
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
