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



}
