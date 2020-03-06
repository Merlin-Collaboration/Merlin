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
#include "MADInterface.h"
#include "ApertureConfiguration.h"

using namespace std;

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

	assert(app->GetType() == "COLLIMATOR");
	assert(app->GetFullEntranceHeight() == 0.1);
	assert(app->GetFullEntranceWidth() == 0.1);
	assert(app->GetCollimatorTilt() == 0);
	assert(app->GetCollimatorLength() == 10);
	assert(app->CheckWithinApertureBoundaries(0.0, 0.0, 0.0) == true);

	collimator_aperture_tilt = 0.001;
	bool side = true;

	OneSidedUnalignedCollimatorAperture* appp = new OneSidedUnalignedCollimatorAperture(x_size_entrance,
		y_size_entrance, collimator_aperture_tilt, length, x_pos_entrance, y_pos_entrance, side);

	assert(appp->GetType() == "ONE-SIDED COLLIMATOR");
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
	DataTable dt;
	dt.AddColumn("APERTYPE", 's');
	dt.AddColumn("S", 'd');
	dt.AddColumn("L", 'd');
	dt.AddColumn("APER_1", 'd');
	dt.AddColumn("APER_2", 'd');
	dt.AddColumn("APER_3", 'd');
	dt.AddColumn("APER_4", 'd');
	dt.AddRow();

	dt.Set("APERTYPE", 0, "RECTELLIPSE");
	dt.Set_d("S", 0, 1);
	dt.Set_d("L", 0, 1);
	dt.Set_d("APER_1", 0, 1);
	dt.Set_d("APER_2", 0, 1);
	dt.Set_d("APER_3", 0, 1);
	dt.Set_d("APER_4", 0, 1);

	DataTableRowIterator itr = dt.begin();
	Aperture* ap = factory.GetInstance(*itr);

	assert(ap->GetType() == "RECTELLIPSE");
	assert(ap->CheckWithinApertureBoundaries(0.0, 0.0, 0.0) == true);
}

void testInterpolatedApertureFactory()
{
	InterpolatorFactory intfactory;
	DataTable dt;
	DataTable dt2;

	dt.AddColumn("APERTYPE", 's');
	dt.AddColumn("S", 'd');
	dt.AddColumn("L", 'd');
	dt.AddColumn("APER_1", 'd');
	dt.AddColumn("APER_2", 'd');
	dt.AddColumn("APER_3", 'd');
	dt.AddColumn("APER_4", 'd');
	dt.AddRow();
	dt.AddRow();

	dt2.AddColumn("APERTYPE", 's');
	dt2.AddColumn("S", 'd');
	dt2.AddColumn("L", 'd');
	dt2.AddColumn("APER_1", 'd');
	dt2.AddColumn("APER_2", 'd');
	dt2.AddColumn("APER_3", 'd');
	dt2.AddColumn("APER_4", 'd');
	dt2.AddRow();
	dt2.AddRow();

	dt.Set("APERTYPE", 0, "RECTELLIPSE");
	dt.Set_d("S", 0, 1);
	dt.Set_d("L", 0, 1);
	dt.Set_d("APER_1", 0, 1);
	dt.Set_d("APER_2", 0, 1);
	dt.Set_d("APER_3", 0, 1);
	dt.Set_d("APER_4", 0, 1);
	dt.Set("APERTYPE", 1, "RECTELLIPSE");
	dt.Set_d("S", 1, 2);
	dt.Set_d("L", 1, 2);
	dt.Set_d("APER_1", 1, 2);
	dt.Set_d("APER_2", 1, 2);
	dt.Set_d("APER_3", 1, 2);
	dt.Set_d("APER_4", 1, 2);

	dt2.Set("APERTYPE", 0, "CIRCLE");
	dt2.Set_d("S", 0, 1);
	dt2.Set_d("L", 0, 1);
	dt2.Set_d("APER_1", 0, 1);
	dt2.Set("APERTYPE", 1, "CIRCLE");
	dt2.Set_d("S", 1, 2);
	dt2.Set_d("L", 1, 2);
	dt2.Set_d("APER_1", 1, 2);

	Aperture* apInt = intfactory.GetInstance(dt);
	Aperture* apInt2 = intfactory.GetInstance(dt2);

	assert(apInt->GetType() == "RECTELLIPSE-interpolated");
	assert(apInt2->GetType() == "CIRCLE-interpolated");
}

int main(int argc, char* argv[])
{
	testApertureFactory();
	testInterpolatedApertureFactory();
	testCollimatorAperture();
	cout << "all aperture tests successful" << endl;
}
