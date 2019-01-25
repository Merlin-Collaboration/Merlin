/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>

#include "InterpolatedApertures.h"
#include "Aperture.h"
#include "CollimatorAperture.h"

void testCollimatorAperture()
{
	double x_size_entrance = 0.1;
	double y_size_entrance = 0.1;
	double collimator_aperture_tilt = 0;
	double length = 10;
	double x_pos_entrance = 0;
	double y_pos_entrance = 0;

	CollimatorAperture* app = new CollimatorAperture(x_size_entrance, y_size_entrance, collimator_aperture_tilt, length,
		x_pos_entrance, y_pos_entrance);

	assert(app->getApertureType() == "COLLIMATOR");
	assert(app->GetFullEntranceHeight() == 0.1);
	assert(app->GetFullEntranceWidth() == 0.1);
	assert(app->GetCollimatorTilt() == 0);
	assert(app->GetCollimatorLength() == 10);
	assert(app->CheckWithinApertureBoundaries(0.0, 0.0, 0.0) == true);

	collimator_aperture_tilt = 0.001;
	bool side = true;

	OneSidedUnalignedCollimatorAperture* appp = new OneSidedUnalignedCollimatorAperture(x_size_entrance,
		y_size_entrance, collimator_aperture_tilt, length, x_pos_entrance, y_pos_entrance, side);

	assert(appp->getApertureType() == "ONE-SIDED COLLIMATOR");
	assert(appp->GetFullEntranceHeight() == 0.1);
	assert(appp->GetFullEntranceWidth() == 0.1);
	assert(appp->GetCollimatorTilt() == 0.001);
	assert(appp->GetCollimatorLength() == 10);
	assert(appp->CheckWithinApertureBoundaries(0.0, 0.0, 0.0) == true);
	assert(appp->GetJawSide() == 1);
	appp->SetJawSide(false);
	assert(appp->GetJawSide() == 0);
}

void testApertureFactory()
{
	ApertureFactory factory;

	string type  = "RECTELLIPSE";
	double s = 1;
	double ap1 = 1;
	double ap2 = 1;
	double ap3 = 1;
	double ap4 = 1;

	Aperture* ap = factory.getInstance(type, s, ap1, ap2, ap3, ap4);

	assert(ap->getRectHalfWidth() == 1);
	assert(ap->getType() == "RECTELLIPSE");
	assert(ap->getApertureType() == "RECTELLIPSE");
	assert(ap->CheckWithinApertureBoundaries(0.0, 0.0, 0.0) == true);
}

void testInterpolatedApertureFactory()
{
	ApertureFactory factory;
	InterpolatorFactory intfactory;

	string type  = "RECTELLIPSE";
	double s = 1;
	double ap1 = 1;
	double ap2 = 1;
	double ap3 = 1;
	double ap4 = 1;

	string type2  = "CIRCLE";
	double s2 = 2;
	double ap12 = 2;
	double ap22 = 2;
	double ap32 = 2;
	double ap42 = 2;

	Aperture* ap = factory.getInstance(type, s, ap1, ap2, ap3, ap4);
	Aperture* app = factory.getInstance(type2, s2, ap12, ap22, ap32, ap42);

	vector<Aperture*> apVec {ap, app};

	Aperture* apInt = intfactory.getInstance(apVec);

	assert(apInt->getType() == "RECTELLIPSEinterpolated");
}

int main(int argc, char* argv[])
{
	testApertureFactory();
	testInterpolatedApertureFactory();
	testCollimatorAperture();
	cout << "all aperture tests successful" << endl;
}
