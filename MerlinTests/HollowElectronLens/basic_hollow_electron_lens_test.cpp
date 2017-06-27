#include "../tests.h"
#include <iostream>
#include <fstream>

#include "AcceleratorModelConstructor.h"
#include "Components.h"
#include "PhysicalUnits.h"
#include "ParticleBunchTypes.h"
#include "ParticleTracker.h"
#include "HollowElectronLens.h"
#include "HollowELensProcess.h"

/*
 * Test that the deterministic modes of the the Hollow Election Lens
 * agree with a reference file.
 *
*/

using namespace std;
using namespace PhysicalUnits;

int main()
{
	AcceleratorModelConstructor* ctor = new AcceleratorModelConstructor();
	ctor->NewModel();
	AcceleratorComponent *hel = new HollowElectronLens("hel1",0);
	ctor->AppendComponent(*hel);
	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;

	const double beam_energy = 7000.0;
	const size_t npart = 400;

	// Test bunch with particles along x, y and diagonals
	vector<Particle> pcoords;
	for(size_t i=0; i<npart; i++)
	{
		double pos = (double(i)-(npart/2)) * 0.1 * millimeter;
		Particle p1(0);
		p1.x() = pos;
		pcoords.push_back(p1);

		Particle p2(0);
		p2.y() = pos;
		pcoords.push_back(p2);

		Particle p3(0);
		p3.x() = pos;
		p3.y() = pos;
		pcoords.push_back(p3);

		Particle p4(0);
		p4.x() = pos;
		p4.y() = -pos;
		pcoords.push_back(p4);
	}

	const string fname = "basic_hollow_electron_lens_test_out.dat";
	ofstream of(fname);
	if (!of)
	{
		cerr << "Could not open: " << fname << endl;
		exit(1);
	}
	of << "# case id x xp y yp" << endl;

	HollowELensProcess* myHELProcess;
	for(int testcase = 0; testcase < 6; testcase++)
	{
		vector<Particle> bcoords{pcoords};
		ProtonBunch* myBunch = new ProtonBunch(beam_energy,1, bcoords);

		int nturns = 1;
		AcceleratorModel::RingIterator ring = theModel->GetRing();
		ParticleTracker* tracker = new ParticleTracker(ring,myBunch);

		cout << "Running case " << testcase << endl;
		switch (testcase)
		{
		case (0):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.195, 2.334948339E4, 3.0);// LHC: 3m, 10KeV, 5A

			// 1 = opposite to protons (focussing)
			myHELProcess->SetElectronDirection(1);
			myHELProcess->SetOpMode(DC);
			myHELProcess->SetRadii(2*millimeter, 5*millimeter);
			break;
		}
		case (1):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.5, 1.5E4, 2.0);// LHC: 3m, 10KeV, 5A

			myHELProcess->SetElectronDirection(1);
			// radial (measured) or perfect profile
			myHELProcess->SetRadialProfile();
			// for LHC hardware we need to scale the radial profile
			myHELProcess->SetLHCRadialProfile();
			myHELProcess->SetOpMode(DC);
			myHELProcess->SetRadii(2*millimeter, 5*millimeter);
			break;
		}
		case (2):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.195, 2.334948339E4, 3.0);// LHC: 3m, 10KeV, 5A

			myHELProcess->SetElectronDirection(0);
			myHELProcess->SetOpMode(DC);
			myHELProcess->SetRadii(2*millimeter, 5*millimeter);
			break;
		}
		case (3):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.7, 2.334948339E4, 3.0);// LHC: 3m, 10KeV, 5A

			myHELProcess->SetElectronDirection(1);
			myHELProcess->SetAC(0.31, .002, 5E-5, 1E3, 2.);
			myHELProcess->SetOpMode(AC);
			myHELProcess->SetRadii(3*millimeter, 5*millimeter);
			break;
		}
		case (4):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.195, 2.0E4, 3.0);// LHC: 3m, 10KeV, 5A

			myHELProcess->SetElectronDirection(0);
			myHELProcess->SetAC(0.31, .002, 5E-5, 1E3, 2.);
			myHELProcess->SetOpMode(AC);
			myHELProcess->SetRadii(3*millimeter, 5*millimeter);
			nturns = 5;
			break;
		}
		case (5):
		{
			myHELProcess = new HollowELensProcess(3, 0, 5, 0.195, 2.334948339E4, 2.0);// LHC: 3m, 10KeV, 5A

			myHELProcess->SetTurnskip(2);
			myHELProcess->SetOpMode(Turnskip);
			myHELProcess->SetRadii(1*millimeter, 4*millimeter);
			nturns = 5;
			break;
		}
		default:
		{
			cerr << "Unknown test case" <<endl;
			abort();
		}
		}
		tracker->AddProcess(myHELProcess);

		for(int n=0; n < nturns; n++)
		{
			cout << "turn: " << n << endl;
			tracker->Track(myBunch);
		}
		size_t i = 0;
		for (auto ip=myBunch->begin(); ip!=myBunch->end(); ip++, i++)
		{
			of << testcase << " " << i << " " <<   ip->x() << " " << ip->xp() << " " <<   ip->y() << " " << ip->yp() <<endl;
		}

		delete tracker;
		delete myBunch;
	}
	of.close();
	return 0;
}
