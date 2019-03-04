/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "AcceleratorModel.h"
#include "SMPTracker.h"
#include "SMPBunch.h"
#include "SMPBunchConstructor.h"
#include "SMPWakeFieldProcess.h"
#include "NormalTransform.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"

#include <fstream>

using namespace SMPTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

using namespace std;

namespace
{

double DispersionCorrectedEmittance(const PSmoments& S)
{
	double s36 = S(ps_Y, ps_DP);
	double s46 = S(ps_YP, ps_DP);
	double dp2 = S.var(ps_DP);
	double s33 = S.var(ps_Y) - s36 * s36 / dp2;
	double s34 = S(ps_Y, ps_YP) - s36 * s46 / dp2;
	double s44 = S.var(ps_YP) - s46 * s46 / dp2;
	return sqrt(s33 * s44 - s34 * s34);
}

}

void OptimumDispersionMatch(AcceleratorModel* accmdl, BeamData beam0)
{
	SMPTracker tracker(accmdl->GetBeamline());
	tracker.AddProcess(new WakeFieldProcess(1));

	ofstream ofs("dispersion_match.dat");
//	beam0.Dyp = 0;

	beam0.Dy = 0.0007;
	beam0.Dyp = -0.00001;
	/***
	    for(double eta_y = 0; eta_y<=0.0008; eta_y +=0.0001)
	        for(double eta_yp = -10.0e-06; eta_yp<=0; eta_yp+=1.0e-06){

	            beam0.Dy = eta_y;
	            beam0.Dyp = eta_yp;
	            cout<<beam0.Dy<<"  "<<beam0.Dyp<<"  "<<flush;
	            ofs<<beam0.Dy<<"  "<<beam0.Dyp<<"  "<<flush;
	 **/
	SMPBunch* bunch = SMPBunchConstructor(beam0, 31, 11).ConstructSMPBunch();

	tracker.Track(bunch);

	PSmoments S;
	bunch->GetMoments(S);
	double p0 = bunch->GetReferenceMomentum();
	p0 *= 1 + S.mean(ps_DP);

	double gamma = p0 / MeV / ElectronMassMeV;
	double geyc = gamma * DispersionCorrectedEmittance(S);

	cout << geyc / nanometer << endl;
	ofs << geyc / nanometer << endl;

	delete bunch;
//	}
}
