#include <iostream>
#include "ScatteringModelsMerlin.h"
#include "Collimators/ScatteringProcess.h"


namespace Collimation
{

ScatteringModelFixed::ScatteringModelFixed(): is_fixed(false)
{
}

void ScatteringModelFixed::AddProcess(Collimation::ScatteringProcess* S)
{
	if (!is_fixed)
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
	for(auto& sp: Processes)
	{
		delete sp;
	}
}


ScatteringModelMerlin::ScatteringModelMerlin()
{
	AddProcess(new ElasticpN());
	AddProcess(new Elasticpn());
	AddProcess(new SingleDiffractive());
	AddProcess(new Rutherford());
	AddProcess(new Inelastic());
	SetScatterType(4); // FIXME, still needed to select EnergyLoss() in CollimateProtonProcess::DoScatter()
	is_fixed = true;
}

ScatteringModelSixTrack::ScatteringModelSixTrack()
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	SetScatterType(0); // FIXME, still needed to select EnergyLoss() in CollimateProtonProcess::DoScatter()
	is_fixed = true;
}

ScatteringModelSixTrackIoniz::ScatteringModelSixTrackIoniz()
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	SetScatterType(1); // FIXME, still needed to select EnergyLoss() in CollimateProtonProcess::DoScatter()
	is_fixed = true;
}

ScatteringModelSixTrackElastic::ScatteringModelSixTrackElastic()
{
	AddProcess(new ElasticpN());
	AddProcess(new Elasticpn());
	AddProcess(new SixTrackSingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	SetScatterType(2); // FIXME, still needed to select EnergyLoss() in CollimateProtonProcess::DoScatter()
	is_fixed = true;
}

ScatteringModelSixTrackSD::ScatteringModelSixTrackSD()
{
	AddProcess(new SixTrackElasticpN());
	AddProcess(new SixTrackElasticpn());
	AddProcess(new SingleDiffractive());
	AddProcess(new SixTrackRutherford());
	AddProcess(new Inelastic());
	SetScatterType(3); // FIXME, still needed to select EnergyLoss() in CollimateProtonProcess::DoScatter()
	is_fixed = true;
}

}
