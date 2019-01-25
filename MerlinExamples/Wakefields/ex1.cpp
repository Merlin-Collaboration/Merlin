/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <fstream>
#include <iostream>

#include "Aperture.h"
#include "AcceleratorModelConstructor.h"
#include "Collimator.h"

#include "BeamData.h"
#include "ParticleBunch.h"
#include "ParticleDistributionGenerator.h"
#include "ParticleTracker.h"

#include "CollimatorWakeProcess.h"
#include "CollimateParticleProcess.h"
#include "CollimatorWakePotentials.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

#include "RandomNG.h"

//    example on using collimator wakefields

// these classes implement examples of wake potentials
// see A.M. Toader et al., EPAC08, Genua, WEPP161 for
// a collection of collimator wakefield formulae
#include "CollimatorPotentialModels.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

int main()
{

	RandomNG::init();

	// initial offset in meters
	double offset = 0.0015;
	// number of modes
	int modes = 5;

//-----------------------------------------------------
//             define the beam
//             SLAC test with 1.19 GeV
//-----------------------------------------------------

	BeamData mybeam;
	mybeam.charge = 2.0e10;
	mybeam.beta_x = 3. * meter;
	mybeam.beta_y = 10. * meter;
	mybeam.emit_x = 0.36 * millimeter;
	mybeam.emit_y = 0.16 * millimeter;
	mybeam.sig_z  = 0.65 * millimeter;
	mybeam.p0     = 1.19 * GeV;

	mybeam.y0 = offset;

	int npart = 10000;

//-----------------------------------------------------
//             Construct bunch
//-----------------------------------------------------
	ParticleBunch* startBunch = new ParticleBunch(npart, NormalParticleDistributionGenerator(), mybeam);
	startBunch->SetCentroid();

//-----------------------------------------------------
// the initial mean and sigma values
//-----------------------------------------------------

	double xmeani = startBunch->GetMoments(0).first;
	double xsigi  = startBunch->GetMoments(0).second;
	double ymeani = startBunch->GetMoments(2).first;
	double ysigi  = startBunch->GetMoments(2).second;

	cout << "Initial bunch parameters:" << endl;
	cout << "(beta_x= " << mybeam.beta_x << ", emit_x=" << mybeam.emit_x << " => sig_x="
		 << sqrt(mybeam.beta_x * mybeam.emit_x) << " [m]" << endl;
	cout << " beta_y=" << mybeam.beta_y << ", emit_y=" << mybeam.emit_y << " => sig_y="
		 << sqrt(mybeam.beta_y * mybeam.emit_y) << " [m]" << endl;
	cout << " offset= " << offset << " [m])" << endl << endl;
	cout << "Mean x =" << xmeani << ' ' << "Sigma x =" << xsigi << endl;
	cout << "Mean y =" << ymeani << ' ' << "Sigma y =" << ysigi << endl;
	cout << "yp angle :" << startBunch->GetMoments(3).first << endl << endl;

//-----------------------------------------------------
//             construct the accelerator model
//             a collimator between 2 drifts
//-----------------------------------------------------

	AcceleratorModelConstructor* accelerator_model = new AcceleratorModelConstructor();
	accelerator_model->NewModel();

	double driftlength1 = 1.0 * meter;
	Drift* drift1       = new Drift("aDrift1", driftlength1);

	double collimatorlength  = 177. * millimeter;
	double collimatorthick   = 0.004 * meter;
	Collimator* collimator      = new Collimator("aCollimator", collimatorlength, collimatorthick);
	double aperturewidth  = 1.9 * millimeter;
	double apertureheight = 1.9 * millimeter;
	Aperture* aperture = new RectangularAperture(aperturewidth, apertureheight);
	collimator->SetAperture(aperture);

	double driftlength2 = 1.0 * meter;
	Drift* drift2       = new Drift("aDrift2", driftlength2);

	accelerator_model->AppendComponent(drift1);
	accelerator_model->AppendComponent(collimator);
	accelerator_model->AppendComponent(drift2);

	// the collimator needs a CollimateParticleProcess
	// and we want to have scattering
	ofstream lossummary("loss_summary.dat");
	CollimateParticleProcess* collimation = new CollimateParticleProcess(0, COLL_AT_EXIT, &lossummary);
	collimation->ScatterAtCollimator(true);

	AcceleratorModel* model = accelerator_model->GetModel();

//-----------------------------------------------------
//             Create a wake potential and
//             add it to the collimator element
//-----------------------------------------------------

// apply the geometric wakefields
	TaperedCollimatorPotentials* collWake
		=  new TaperedCollimatorPotentials(modes, aperturewidth / 2, apertureheight / 2);
	collimator->SetWakePotentials(collWake);

// //another example: a resistive wakefield
// double conductivity = 2.38e6;         //titanium
// double conductivity = 5.98*10000000;  //copper
// double conductivity = 3.08e7;         //berrilium
// double conductivity = 6.e4;           //carbon
// double conductivity = 4.5e6;          //for TiN
// ResistiveWakePotentials* resWake =  new ResistiveWakePotentials(modes, aperturewidth/2, conductivity, collimatorlength);
// collimator->SetWakePotentials(resWake);

//-----------------------------------------------------
//             Create the WakeProcess
//             and tie it to the wake potential
//-----------------------------------------------------

	CollimatorWakeProcess* collWakeProc = new CollimatorWakeProcess(modes, 1, 100, 3);

	// here we connect the CollimatorWakeProcess
	// and the TaperedCollimatorPotentials
	// i.e. the CollimatorWakeProcess will only
	// take care if he finds a wake potential
	// of type TaperedCollimatorPotentials
	collWake->SetExpectedProcess(collWakeProc);

//-----------------------------------------------------
//             Construct Tracker and add
//             the different processes
//-----------------------------------------------------

	ParticleTracker* tracker = new ParticleTracker(model->GetBeamline(), startBunch);

	tracker->AddProcess(collimation);
	tracker->AddProcess(collWakeProc);

//-----------------------------------------------------
//             Final bunch parameters
//-----------------------------------------------------

	ParticleBunch* finalBunch = tracker->Track(startBunch);

	ofstream final ("FinalParameters.dat");
	finalBunch->Output(final);

	int i = 0;
	double averageyp = 0;
	for(ParticleBunch::iterator p = finalBunch->begin(); p != finalBunch->end(); p++)
	{
		averageyp += (p->yp() - averageyp) / (++i);
	}

	double xmeanf = finalBunch->GetMoments(0).first;
	double xsigf  = finalBunch->GetMoments(0).second;
	double ymeanf = finalBunch->GetMoments(2).first;
	double ysigf  = finalBunch->GetMoments(2).second;

	cout << endl << "Final bunch parameters:" << endl;
	cout << "Mean x =" << xmeanf << ' ' << "Sigma x =" << xsigf << endl;
	cout << "Mean y =" << ymeanf << ' ' << "Sigma y =" << ysigf << endl;
	cout << "yp angle final:" << finalBunch->GetMoments(3).first << endl;

	return 0;

}
