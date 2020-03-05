/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include "ScatteringModelsMerlin.h"
#include "ScatteringProcess.h"
#include "ElasticScatter.h"
#include "DiffractiveScatter.h"

using namespace PhysicalConstants;

namespace Collimation
{

void ScatteringModelMerlin::Configure(MaterialProperties * m, double Energy)
{
	Processes[1] = new Rutherford(m);
	Processes[2] = new Elasticpn(Energy);
	Processes[3] = new SingleDiffractive(Energy);
	Processes[4] = new Inelastic();
	Processes[5] = new ElasticpN(Energy, m);
	Xsection[0] = m->sigma_T;
	Xsection[1] = m->sigma_R;
	Xsection[2] = 1.618 * pow(m->A, 0.333) * Processes[2]->sigma;
	Xsection[3] = 1.618 * pow(m->A, 0.333) * Processes[3]->sigma;
	Xsection[4] = m->sigma_I;
	std::cout << "'Merlin' Cross sections Total " << Xsection[0] << " Rutherford " << Xsection[1] << " Elastic  "
			  << Xsection[2] << " Diffractive "
			  << Xsection[3] << " Inelastic  " << Xsection[4] << std::endl;
}

void ScatteringModelSixTrack::Configure(MaterialProperties * m, double Energy)
{
	double s = 2 * ProtonMassGeV * Energy + ProtonMassGeV * ProtonMassGeV;
	Processes[1] = new SixTrackRutherford();
	Processes[2] = new SixTrackElasticpn();
	Processes[3] = new SixTrackSingleDiffractive();
	Processes[4] = new Inelastic();
	Processes[5] = new SixTrackElasticpN(m);
	Xsection[0] = m->sigma_T;
	Xsection[1] = m->sigma_R;
	Xsection[2] = 1.618 * pow(m->A, 0.333) * 0.007 * pow(Energy / 450.0, 0.04792);
	Xsection[3] = 1.618 * pow(m->A, 0.333) * 0.00068 * log(0.15 * s);
	Xsection[4] = m->sigma_I;
	energy_loss_mode = SimpleEnergyLoss;
	std::cout << "'Sixtrack' cross sections Total " << Xsection[0] << " Rutherford " << Xsection[1] << " Elastic  "
			  << Xsection[2] << " Diffractive "
			  << Xsection[3] << " Inelastic  " << Xsection[4] << std::endl;
}

void ScatteringModelSixTrackIoniz::Configure(MaterialProperties * m, double Energy)
{
	double s = 2 * ProtonMassGeV * Energy + ProtonMassGeV * ProtonMassGeV;
	Processes[1] = new SixTrackRutherford();
	Processes[2] = new SixTrackElasticpn();
	Processes[3] = new SixTrackSingleDiffractive();
	Processes[4] = new Inelastic();
	Processes[5] = new SixTrackElasticpN(m);
	Xsection[0] = m->sigma_T;
	Xsection[1] = m->sigma_R;
	Xsection[2] = 1.618 * pow(m->A, 0.333) * 0.007 * pow(Energy / 450.0, 0.04792);
	Xsection[3] = 1.618 * pow(m->A, 0.333) * 0.00068 * log(0.15 * s);
	Xsection[4] = m->sigma_I;
	energy_loss_mode = FullEnergyLoss;
	std::cout << "'Sixtrack' cross sections Total " << Xsection[0] << " Rutherford " << Xsection[1] << " Elastic  "
			  << Xsection[2] << " Diffractive "
			  << Xsection[3] << " Inelastic  " << Xsection[4] << std::endl;
}

void ScatteringModelSixTrackElastic::Configure(MaterialProperties * m, double Energy)
{
	double s = 2 * ProtonMassGeV * Energy + ProtonMassGeV * ProtonMassGeV;
	Processes[1] = new SixTrackRutherford();
	Processes[2] = new Elasticpn(Energy);
	Processes[3] = new SixTrackSingleDiffractive();
	Processes[4] = new Inelastic();
	Processes[5] = new ElasticpN(Energy, m);
	Xsection[0] = m->sigma_T;
	Xsection[1] = m->sigma_R;
	Xsection[2] = 1.618 * pow(m->A, 0.333) * Processes[2]->sigma;
	Xsection[3] = 1.618 * pow(m->A, 0.333) * 0.00068 * log(0.15 * s);
	Xsection[4] = m->sigma_I;
	energy_loss_mode = SimpleEnergyLoss;
	std::cout << "'Sixtrack' cross sections Total " << Xsection[0] << " Rutherford " << Xsection[1] << " Elastic  "
			  << Xsection[2] << " Diffractive "
			  << Xsection[3] << " Inelastic  " << Xsection[4] << std::endl;
}

void ScatteringModelSixTrackSD::Configure(MaterialProperties * m, double Energy)
{
	Processes[1] = new SixTrackRutherford();
	Processes[2] = new SixTrackElasticpn();
	Processes[3] = new SingleDiffractive(Energy);
	Processes[4] = new Inelastic();
	Processes[5] = new SixTrackElasticpN(m);
	Xsection[0] = m->sigma_T;
	Xsection[1] = m->sigma_R;
	Xsection[2] = 1.618 * pow(m->A, 0.333) * 0.007 * pow(Energy / 450.0, 0.04792);
	Xsection[3] = 1.618 * pow(m->A, 0.333) * Processes[3]->sigma;
	Xsection[4] = m->sigma_I;
	energy_loss_mode = SimpleEnergyLoss;
	std::cout << "'Sixtrack' cross sections Total " << Xsection[0] << " Rutherford " << Xsection[1] << " Elastic  "
			  << Xsection[2] << " Diffractive "
			  << Xsection[3] << " Inelastic  " << Xsection[4] << std::endl;
}

}
