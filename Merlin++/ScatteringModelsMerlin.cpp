/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include "ScatteringModelsMerlin.h"
#include "ScatteringProcess.h"

namespace Collimation
{

ScatteringModelFixed::ScatteringModelFixed(MaterialProperties* mat) :
	ScatteringModel(),  is_fixed(false)
{
}

void ScatteringModelFixed::AddProcess(Collimation::ScatteringProcess* S)
{
	if(!is_fixed)
	{
		ScatteringModel::AddProcess(S);
	}
	else
	{
		std::cerr << "Can't add mode processes to a fixed ScatteringModel" << std::endl;
		abort();
	}
}

ScatteringModelFixed::~ScatteringModelFixed()
{
	// Processes can only be added by the class itself (as long as the constructor
	// sets is_fixed=true) so they can be safely deleted by the class.
	for(auto& sp : Processes)
	{
		delete sp;
	}
}

ScatteringModelMerlin::ScatteringModelMerlin(MaterialProperties* mat) :
	ScatteringModelFixed(mat)
{
	AddProcess(new ElasticpN(7000));
	AddProcess(new Elasticpn(7000.0));
	AddProcess(new SingleDiffractive(7000.0));
	AddProcess(new Rutherford(mat));
	AddProcess(new Inelastic());
	SetScatterType(4); // FIXME, still needed for CrossSections
	is_fixed = true;
}

ScatteringModelSixTrack::ScatteringModelSixTrack(MaterialProperties* mat) :
	ScatteringModelFixed(mat)
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	energy_loss_mode = SimpleEnergyLoss;
	SetScatterType(0); // FIXME, still needed for CrossSections
	is_fixed = true;
}

ScatteringModelSixTrackIoniz::ScatteringModelSixTrackIoniz(MaterialProperties* mat) :
	ScatteringModelFixed(mat)
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	SetScatterType(1); // FIXME, still needed for CrossSections
	is_fixed = true;
}

ScatteringModelSixTrackElastic::ScatteringModelSixTrackElastic(MaterialProperties*  mat) :
	ScatteringModelFixed(mat)
{
	AddProcess(new ElasticpN(7000.0));
	AddProcess(new Elasticpn(7000.0));
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	energy_loss_mode = SimpleEnergyLoss;
	SetScatterType(2); // FIXME, still needed for CrossSections
	is_fixed = true;
}

ScatteringModelSixTrackSD::ScatteringModelSixTrackSD(MaterialProperties* mat) :
	ScatteringModelFixed(mat)
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SingleDiffractive(7000.0));
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	energy_loss_mode = SimpleEnergyLoss;
	SetScatterType(3); // FIXME, still needed for CrossSections
	is_fixed = true;
}

}
