/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include <iostream>
#include "../tests.h"
#include "MaterialData.h"
#include "PhysicalUnits.h"
#include "RandomNG.h"

using namespace std;
using namespace PhysicalUnits;

int main(int argc, char* argv[])
{
	RandomNG::init();
	auto md = StandardMaterialData();

	cout << "Al Z  : " << md.property["Al"]->Z << endl;
	assert(md.property["Al"]->Z == 13);
	cout << "Al MEE: " << md.property["Al"]->HaveExtra("MeanExcitationEnergy")
		 << " " << md.property["Al"]->GetExtra("MeanExcitationEnergy") << endl;
	assert(md.property["Al"]->HaveExtra("MeanExcitationEnergy"));
	cout << "Al    : " << *(md.property["Al"]) << endl;

	md.MakeMixture("rust", "Fe O", 2.0, 3.0, 5.25 * gram / cc);
	md.MakeMixtureByWeight("rust2", "Fe O", 111.69, 47.9982, 5.25 * gram / cc);

	md.MakeMixture("rust3", "Fe O", {2.0, 3.0}, 5.25 * gram / cc);
	md.MakeMixtureByWeight("rust4", "Fe O", {111.69, 47.9982}, 5.25 * gram / cc);

	vector<double> props;
	props.push_back(2);
	props.push_back(3);
	md.MakeMixture("rust5", "Fe O", props, 5.25 * gram / cc);

	md.property["rust"]->SetExtra("HeatCap", 103.9);

	cout << "rust  : " << *(md.property["rust"]) << endl;
	cout << "rust2 : " << *(md.property["rust2"]) << endl;
	cout << "rust3 : " << *(md.property["rust3"]) << endl;
	cout << "rust4 : " << *(md.property["rust4"]) << endl;

	assert_close(md.property["rust"]->A, 31.93764, 1e-3);
	assert_close(md.property["rust"]->Z, 15.2, 1e-3);
	assert_close(md.property["rust2"]->A, 31.93764, 1e-3);
	assert_close(md.property["rust3"]->A, 31.93764, 1e-3);
	assert_close(md.property["rust4"]->A, 31.93764, 1e-3);
	assert_close(md.property["rust5"]->A, 31.93764, 1e-3);

	md.PrintTable();

	// check that A_H() and A_R() give value for one of the components
	map<double, int> count_H, count_R;
	for(int i = 0; i < 100; i++)
	{
		count_H[md.property["rust"]->A_H()] += 1;
		count_R[md.property["rust"]->A_R()] += 1;
	}

	assert(count_H[md.property["O"]->A] > 0);
	assert(count_H[md.property["Fe"]->A] > 0);
	assert(count_H[md.property["O"]->A] + count_H[md.property["Fe"]->A] == 100);

	//assert(count_R[md.property["O"]->A] > 0); // O sigma_T is -1
	assert(count_R[md.property["Fe"]->A] > 0);
	assert(count_R[md.property["O"]->A] + count_R[md.property["Fe"]->A] == 100);

	return 0;
}
