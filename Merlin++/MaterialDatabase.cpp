/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>

#include "MaterialDatabase.h"
#include "MaterialMixture.h"
#include "CompositeMaterial.h"

#include "PhysicalUnits.h"

using namespace PhysicalUnits;

MaterialDatabase::MaterialDatabase()
{
	/*
	   Here we create new materials, add their properties, then push them into a map for manipulation.
	   References: Particle data group: http://pdg.lbl.gov/2013/AtomicNuclearProperties/
	 */

//ALL CROSS SECTIONS IN BARNS!

	//new Material(Name, A, AN, SE, SI, SR, dE, X, rho, sig);
	//Beryllium - using constructor rather than addinf each variable
	Material* Be = new Material("Beryllium", "Be", 9.012182, 4, 0.069, 0.199, 0.000035, 0.55, 651900, 1848, 3.08E7);
	Be->SetSixtrackTotalNucleusCrossSection(0.268);
	Be->SetSixtrackNuclearSlope(74.7);
	Be->SetMeanExcitationEnergy(63.7 * eV);
	Be->SetElectronDensity(Be->CalculateElectronDensity());
	Be->SetPlasmaEnergy(Be->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Be->GetSymbol(), Be));

	//Carbon (graphite)
	Material* C = new Material();
	C->SetAtomicNumber(6);
	C->SetAtomicMass(12.0107);
	C->SetName("Carbon");
	C->SetSymbol("C");
	C->SetSixtrackTotalNucleusCrossSection(0.331);
	C->SetSixtrackInelasticNucleusCrossSection(0.231);
	C->SetSixtrackRutherfordCrossSection(0.000076);
	C->SetSixtrackdEdx(0.68);
	C->SetConductivity(7.14E4);
	C->SetRadiationLength(427000);
	C->SetDensity(2210);
	C->SetSixtrackNuclearSlope(70.0);
	C->SetMeanExcitationEnergy(78.0 * eV);
	C->SetElectronDensity(C->CalculateElectronDensity());
	C->SetPlasmaEnergy(C->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(C->GetSymbol(), C));

	//Oxygen (note the density, etc is for liquid)!
	//Some of these parameters are rather invalid
	//Also note that bits of the scattering code currently are invalid for non-solids.
	Material* O = new Material();
	O->SetAtomicNumber(8);
	O->SetAtomicMass(15.9994);
	O->SetName("Oxygen");
	O->SetSymbol("O");
	O->SetSixtrackTotalNucleusCrossSection(O->CalculateSixtrackTotalNucleusCrossSection());
	O->SetSixtrackInelasticNucleusCrossSection(O->CalculateSixtrackInelasticNucleusCrossSection());
	O->SetSixtrackRutherfordCrossSection(0.0133766);
	O->SetSixtrackdEdx(O->CalculateSixtrackdEdx());
	O->SetConductivity(1);  //See here
	O->SetRadiationLength(O->CalculateRadiationLength());
	O->SetDensity(1140);    //And here
	O->SetSixtrackNuclearSlope(O->CalculateSixtrackNuclearSlope());
	O->SetMeanExcitationEnergy(95.0 * eV);
	O->SetElectronDensity(O->CalculateElectronDensity());
	O->SetPlasmaEnergy(O->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(O->GetSymbol(), O));

	//Aluminium
	Material* Al = new Material();
	Al->SetAtomicNumber(13);
	Al->SetAtomicMass(26.9815386);
	Al->SetName("Aluminium");
	Al->SetSymbol("Al");
	Al->SetSixtrackTotalNucleusCrossSection(0.634);
	Al->SetSixtrackInelasticNucleusCrossSection(0.421);
	Al->SetSixtrackRutherfordCrossSection(0.00034);
	Al->SetSixtrackdEdx(0.81);
	Al->SetConductivity(3.564E7);
	Al->SetRadiationLength(Al->CalculateRadiationLength());
	Al->SetDensity(2699);
	Al->SetSixtrackNuclearSlope(120.3);
	Al->SetMeanExcitationEnergy(166.0 * eV);
	Al->SetElectronDensity(Al->CalculateElectronDensity());
	Al->SetPlasmaEnergy(Al->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Al->GetSymbol(), Al));

	//Iron
	Material* Fe = new Material();
	Fe->SetAtomicNumber(26);
	Fe->SetAtomicMass(55.845);
	Fe->SetName("Iron");
	Fe->SetSymbol("Fe");
	Fe->SetSixtrackTotalNucleusCrossSection(Fe->CalculateSixtrackTotalNucleusCrossSection());
	Fe->SetSixtrackInelasticNucleusCrossSection(Fe->CalculateSixtrackInelasticNucleusCrossSection());
	Fe->SetSixtrackRutherfordCrossSection(Fe->CalculateSixtrackRutherfordCrossSection());
	Fe->SetSixtrackdEdx(Fe->CalculateSixtrackdEdx());
	Fe->SetConductivity(1.04E7);
	Fe->SetRadiationLength(Fe->CalculateRadiationLength());
	Fe->SetDensity(7870);
	Fe->SetSixtrackNuclearSlope(Fe->CalculateSixtrackNuclearSlope());
	Fe->SetMeanExcitationEnergy(286.0 * eV);
	Fe->SetElectronDensity(Fe->CalculateElectronDensity());
	Fe->SetPlasmaEnergy(Fe->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Fe->GetSymbol(), Fe));

	//Nickel
	Material* Ni = new Material();
	Ni->SetAtomicNumber(28);
	Ni->SetAtomicMass(58.6934);
	Ni->SetName("Nickel");
	Ni->SetSymbol("Ni");
	Ni->SetSixtrackTotalNucleusCrossSection(Ni->CalculateSixtrackTotalNucleusCrossSection());
	Ni->SetSixtrackInelasticNucleusCrossSection(Ni->CalculateSixtrackInelasticNucleusCrossSection());
	Ni->SetSixtrackRutherfordCrossSection(Ni->CalculateSixtrackRutherfordCrossSection());
	Ni->SetSixtrackdEdx(Ni->CalculateSixtrackdEdx());
	Ni->SetConductivity(1.44E7);
	Ni->SetRadiationLength(Ni->CalculateRadiationLength());
	Ni->SetDensity(8900);
	Ni->SetSixtrackNuclearSlope(Ni->CalculateSixtrackNuclearSlope());
	Ni->SetMeanExcitationEnergy(311.0 * eV);
	Ni->SetElectronDensity(Ni->CalculateElectronDensity());
	Ni->SetPlasmaEnergy(Ni->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Ni->GetSymbol(), Ni));

	//Copper
	Material* Cu = new Material();
	Cu->SetAtomicNumber(29);
	Cu->SetAtomicMass(63.546);
	Cu->SetName("Copper");
	Cu->SetSymbol("Cu");
	Cu->SetSixtrackTotalNucleusCrossSection(1.232);
	Cu->SetSixtrackInelasticNucleusCrossSection(0.782);
	Cu->SetSixtrackRutherfordCrossSection(0.00153);
	Cu->SetSixtrackdEdx(2.69);
	Cu->SetConductivity(5.98E7);
	Cu->SetRadiationLength(Cu->CalculateRadiationLength());
	Cu->SetDensity(8960);
	Cu->SetSixtrackNuclearSlope(217.8);
	Cu->SetMeanExcitationEnergy(322.0 * eV);
	Cu->SetElectronDensity(Cu->CalculateElectronDensity());
	Cu->SetPlasmaEnergy(Cu->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Cu->GetSymbol(), Cu));

	//Molybdenum
	Material* Mo = new Material();
	Mo->SetAtomicNumber(42);
	Mo->SetAtomicMass(95.96);
	Mo->SetName("Molybdenum");
	Mo->SetSymbol("Mo");
	Mo->SetSixtrackTotalNucleusCrossSection(Mo->CalculateSixtrackTotalNucleusCrossSection());
	Mo->SetSixtrackInelasticNucleusCrossSection(Mo->CalculateSixtrackInelasticNucleusCrossSection());
	Mo->SetSixtrackRutherfordCrossSection(0.264483);
	Mo->SetSixtrackdEdx(Mo->CalculateSixtrackdEdx());
	Mo->SetConductivity(1.87E7);
	Mo->SetRadiationLength(Mo->CalculateRadiationLength());
	Mo->SetDensity(10200);
	Mo->SetSixtrackNuclearSlope(Mo->CalculateSixtrackNuclearSlope());
	Mo->SetMeanExcitationEnergy(424.0 * eV);
	Mo->SetElectronDensity(Mo->CalculateElectronDensity());
	Mo->SetPlasmaEnergy(Mo->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Mo->GetSymbol(), Mo));

	//Tungsten
	Material* W = new Material();
	W->SetAtomicNumber(74);
	W->SetAtomicMass(183.84);
	W->SetName("Tungsten");
	W->SetSymbol("W");
	W->SetSixtrackTotalNucleusCrossSection(2.767);
	W->SetSixtrackInelasticNucleusCrossSection(1.65);
	W->SetSixtrackRutherfordCrossSection(0.00768);
	W->SetSixtrackdEdx(5.79);
	W->SetConductivity(1.77E3);
	W->SetRadiationLength(W->CalculateRadiationLength());
	W->SetDensity(19300);
	W->SetSixtrackNuclearSlope(440.3);
	W->SetMeanExcitationEnergy(727.0 * eV);
	W->SetElectronDensity(W->CalculateElectronDensity());
	W->SetPlasmaEnergy(W->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(W->GetSymbol(), W));

	//Lead
	Material* Pb = new Material();
	Pb->SetAtomicNumber(82);
	Pb->SetAtomicMass(207.2);
	Pb->SetName("Lead");
	Pb->SetSymbol("Pb");
	Pb->SetSixtrackTotalNucleusCrossSection(2.960);
	Pb->SetSixtrackInelasticNucleusCrossSection(1.77);
	Pb->SetSixtrackRutherfordCrossSection(0.00907);
	Pb->SetSixtrackdEdx(3.40);
	Pb->SetConductivity(4.8077E6);
	Pb->SetRadiationLength(Pb->CalculateRadiationLength());
	Pb->SetDensity(11350);
	Pb->SetSixtrackNuclearSlope(455.3);
	Pb->SetMeanExcitationEnergy(823.0 * eV);
	Pb->SetElectronDensity(Pb->CalculateElectronDensity());
	Pb->SetPlasmaEnergy(Pb->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(Pb->GetSymbol(), Pb));

	/*
	 * Other assorted materials
	 */

	//From H.Day thesis + CERN-ATS-2011-224
	//CuDiamond Conductivity = 1.25E7	density ~5400
	//MoDiamond Conductivity = 5.5E6	density ~6900

	/*
	 * AC150K
	 * LHC TCP/TCSG collimator material
	 * Graphite but different density etc.
	 */
	Material* AC150K = new Material();
	AC150K->SetAtomicNumber(6);
	AC150K->SetAtomicMass(12.0107);
	AC150K->SetName("AC150K");
	AC150K->SetSymbol("AC150K");
	AC150K->SetSixtrackTotalNucleusCrossSection(0.331);
	AC150K->SetSixtrackInelasticNucleusCrossSection(0.231);
	AC150K->SetSixtrackRutherfordCrossSection(0.000076);
	AC150K->SetSixtrackdEdx(0.68);      //Needs to be scaled
	AC150K->SetConductivity(7.14E4);    //This is just the same value as graphite - check
	AC150K->SetRadiationLength(AC150K->CalculateRadiationLength());
	AC150K->SetDensity(1650);       //This is different - was 2210
	AC150K->SetSixtrackNuclearSlope(70.0);
	AC150K->SetMeanExcitationEnergy(78.0 * eV);
	AC150K->SetElectronDensity(AC150K->CalculateElectronDensity());
	AC150K->SetPlasmaEnergy(AC150K->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(AC150K->GetSymbol(), AC150K));

	/*
	 * INERMET 180
	 * LHC TCLA/TCT material
	 * Tungsten alloy with Nickel and Copper
	 */
	MaterialMixture* IT180 = new MaterialMixture();
	IT180->SetName("INERMET180");
	IT180->SetSymbol("IT180");
	IT180->AddMaterialByMassFraction(W, 0.95);
	IT180->AddMaterialByMassFraction(Ni, 0.035);
	IT180->AddMaterialByMassFraction(Cu, 0.015);
	IT180->Assemble();

	IT180->SetDensity(18060);
	IT180->SetConductivity(1.77E3); //Set to the same as Tungsten for now
	IT180->SetSixtrackdEdx(IT180->CalculateSixtrackdEdx());
	IT180->SetRadiationLength(IT180->CalculateRadiationLength());
	IT180->SetMeanExcitationEnergy(IT180->CalculateMeanExcitationEnergy());
	IT180->SetElectronDensity(IT180->CalculateElectronDensity());
	IT180->SetPlasmaEnergy(IT180->CalculatePlasmaEnergy());
	db.insert(pair<string, Material*>(IT180->GetSymbol(), IT180));

	/*
	 * LHC TCLA/TCL material
	 * Copper + Aluminium Oxide powder
	 * Glidcop - 99.72% Copper + 0.28% Aluminium Oxide for Glidcop-15 by mass
	 */
	//~ MaterialMixture* Glidcop = new MaterialMixture();
	CompositeMaterial* Glidcop = new CompositeMaterial();
	Glidcop->SetName("Glidcop");
	Glidcop->SetSymbol("GCOP");
	double Al_M = 0.4 * Al->GetAtomicNumber();
	double O_M = 0.6 * O->GetAtomicNumber();
	Glidcop->AddMaterialByMassFraction(Cu, 0.9972);
	Glidcop->AddMaterialByMassFraction(Al, 0.0028 * Al_M / (Al_M + O_M));
	Glidcop->AddMaterialByMassFraction(O, 0.0028 * O_M / (Al_M + O_M));

	Glidcop->SetDensity(8930);
	Glidcop->SetConductivity(5.38E7);   //CERN-ATS-2011-224

	Glidcop->Assemble();

	// Test Mixture function
	std::vector<std::pair<std::string, double> > els = Glidcop->GetConstituentElements();
	std::vector<std::pair<std::string, double> >::iterator el_it;

	db.insert(std::pair<std::string, Material*>(Glidcop->GetSymbol(), Glidcop));
}

//Try and find the material we want
Material* MaterialDatabase::FindMaterial(std::string symbol)
{
	//Iterator for use with the material pointers map
	std::map<std::string, Material*>::iterator position;

	//Try to find the material required
	position = db.find(symbol);

	//did not find the material, this is bad.
	if(position == db.end())
	{
		std::cerr << "Requested aperture material '" << symbol << "' not found. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	//This should return a pointer to the material we are interested in.
	return position->second;
}

bool MaterialDatabase::VerifyMaterials()
{
	bool verification = true;
	std::map<string, Material*>::const_iterator MaterialIt;
	MaterialIt = db.begin();
	while(MaterialIt != db.end())
	{
		bool test = MaterialIt->second->VerifyMaterial();
		if(test == false)
		{
			verification = false;
			std::cerr << "Material " << MaterialIt->second->GetName() << " failed to verify." << std::endl;
		}
		MaterialIt++;
	}
	return verification;
}

void MaterialDatabase::DumpMaterialProperties()
{
	std::cout.precision(5);
	std::cout << std::setw(15) << "Name" << std::setw(7) << "Symb" << std::setw(4) << "Z";
	std::cout << std::setw(7) << "A" << std::setw(6) << "dnsty" << std::setw(7) << "X0 mm";
	std::cout << std::setw(10) << "X0 m";
	std::cout << std::setw(11) << "n_e" << std::setw(11) << "E_plas";
	std::cout << std::setw(8) << "sig_tot" << std::setw(8) << "sig_in" << std::setw(8) << "sig_R";
	std::cout << std::setw(8) << "dEdx";
	std::cout << std::endl;

	std::map<string, Material*>::const_iterator MaterialIt;
	MaterialIt = db.begin();
	while(MaterialIt != db.end())
	{
		Material* m = MaterialIt->second;
		if(!m->IsMixture())
		{
			std::cout << std::setw(15) << m->GetName();
			std::cout << std::setw(7) << m->GetSymbol();
			std::cout << std::setw(4) << m->GetAtomicNumber();
			std::cout << std::setw(7) << m->GetAtomicMass();
			std::cout << std::setw(6) << m->GetDensity();
			std::cout << std::setw(7) << m->GetRadiationLength() / m->GetDensity();
			std::cout << std::setw(10) << m->GetRadiationLengthInM();
			std::cout << std::setw(11) << m->GetElectronDensity();
			std::cout << std::setw(11) << m->GetPlasmaEnergy();

			std::cout << std::setw(8) << m->GetSixtrackTotalNucleusCrossSection();
			std::cout << std::setw(8) << m->GetSixtrackInelasticNucleusCrossSection();
			std::cout << std::setw(8) << m->GetSixtrackRutherfordCrossSection();
			std::cout << std::setw(8) << m->GetSixtrackdEdx();
			std::cout << std::endl;
		}
		MaterialIt++;
	}

	std::cout << "Mixtures:" <<  std::endl;
	std::cout << std::setw(15) << "Name" << std::setw(7) << "Symb";
	std::cout << std::setw(6) << "dnsty" << std::setw(7) << "X0 mm";
	std::cout << std::setw(10) << "X0 m";
	std::cout << std::setw(11) << "n_e" << std::setw(11) << "E_plas";
	std::cout << std::setw(8) << "dEdx";
	std::cout << std::endl;
	MaterialIt = db.begin();
	while(MaterialIt != db.end())
	{
		Material* m = MaterialIt->second;
		if(m->IsMixture())
		{
			std::cout << std::setw(15) << m->GetName();
			std::cout << std::setw(7) << m->GetSymbol();
			std::cout << std::setw(6) << m->GetDensity();
			std::cout << std::setw(7) << m->GetRadiationLength() / m->GetDensity();
			std::cout << std::setw(10) << m->GetRadiationLengthInM();
			std::cout << std::setw(11) << m->GetElectronDensity();
			std::cout << std::setw(11) << m->GetPlasmaEnergy();
			std::cout << std::setw(8) << m->GetSixtrackdEdx();
			std::cout << std::endl;
		}
		MaterialIt++;
	}
}
