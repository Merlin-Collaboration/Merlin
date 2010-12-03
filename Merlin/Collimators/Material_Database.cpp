#include "Collimators/Material_Database.hpp"
#include "Collimators/Material.hpp"
#include "NumericalUtils/PhysicalUnits.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

using namespace std;
using namespace PhysicalUnits;

Material_Database::Material_Database() : db(0) //size of the db vector, aka number of materials, increase as required.
{

//Here we create new materials, add their properties, then push them onto a vector for manipulation.
//Need to get references for all data.
/*
References:
Particle data group: http://pdg.lbl.gov/2010/AtomicNuclearProperties/
*/

//ALL CROSS SECTIONS IN BARNS

	//Beryllium
	material* Be = new material();
	Be->atomic_number	=	4;
	Be->name		=	"Beryllium";
	Be->symbol		=	"Be";
	Be->sigma_pN_total	=	0.268;
	Be->sigma_pN_inelastic	=	0.199;
	Be->sigma_Rutherford	=	0.000035;
	Be->dEdx		=	0.55;
	Be->rho			=	1.848;
	Be->A			=	9.012182;
	Be->sigma		=	3.08E7;
	Be->X0			=	0.3528;
	Be->density		=	1848;
	Be->b_N			=	74.7;
	Be->MeanExcitationEnergy=	63.7*eV;
	Be->ElectronDensity	=	Be->CalculateElectronDensity();
	Be->PlasmaEnergy	=	Be->CalculatePlasmaEnergy();
	db.push_back(Be);

	//Carbon (graphite)
	material* C = new material();
	C->atomic_number	=	6;
	C->name			=	"Carbon";
	C->symbol		=	"C";
	C->sigma_pN_total	=	0.331;
	C->sigma_pN_inelastic	=	0.231;
	C->sigma_Rutherford	=	0.000076;
	C->dEdx			=	0.68;
	C->rho			=	2.265;		//Check
	C->A			=	12.0107;
	C->sigma		=	7.14E4;
	//C->X0			=	18.8E-2;	//Not for graphite
	C->X0			=	0.1932;
	C->density		=	2210;
	C->b_N			=	70.0;
	C->MeanExcitationEnergy	=	78.0*eV;
	C->ElectronDensity	=	C->CalculateElectronDensity();
	C->PlasmaEnergy		=	C->CalculatePlasmaEnergy();
	db.push_back(C);

	//Aluminium
	material* Al = new material();
	Al->atomic_number	=	13;
	Al->name		=	"Aluminium";
	Al->symbol		=	"Al";
	Al->sigma_pN_total	=	0.634;
	Al->sigma_pN_inelastic	=	0.421;
	Al->sigma_Rutherford	=	0.00034;
	Al->dEdx		=	0.81;
	Al->rho			=	2.70;
	Al->A			=	26.9815386;
	Al->sigma		=	35.64E6;
	Al->X0			=	0.08897;
	Al->density		=	2699;
	Al->b_N			=	120.3;
	Al->MeanExcitationEnergy=	166.0*eV;
	Al->ElectronDensity	=	Al->CalculateElectronDensity();
	Al->PlasmaEnergy	=	Al->CalculatePlasmaEnergy();
	db.push_back(Al);

	//Copper
	material* Cu = new material();
	Cu->atomic_number		=	29;
	Cu->name			=	"Copper";
	Cu->symbol			=	"Cu";
	Cu->sigma_pN_total		=	1.232;
	Cu->sigma_pN_inelastic		=	0.782;
	Cu->sigma_Rutherford		=	0.00153;
//	Cu->dEdx			=	2.69;
	Cu->dEdx			=	1.250776630157339;
	Cu->rho				=	8.96;
	Cu->A				=	63.546;
	Cu->sigma			=	5.98E7;
	Cu->X0				=	0.01436;
	Cu->density			=	8960;
	Cu->b_N				=	217.8;
	Cu->ElectronCriticalEnergy	=	19.42*MeV;
	Cu->MeanExcitationEnergy	=	322.0*eV;
	Cu->ElectronDensity		=	Cu->CalculateElectronDensity();
	Cu->PlasmaEnergy		=	Cu->CalculatePlasmaEnergy();
	db.push_back(Cu);

	//Tungsten
	material* W = new material();
	W->atomic_number	=	74;
	W->name			=	"Tungsten";
	W->symbol		=	"W";
	W->sigma_pN_total	=	2.767;
	W->sigma_pN_inelastic	=	1.65;
	W->sigma_Rutherford	=	0.00768;
	W->dEdx			=	5.79;
	W->rho			=	19.3;
	W->A			=	183.84;
	W->sigma		=	0.177E4;
	W->X0			=	0.003504;
	W->density		=	19300;
	W->b_N			=	440.3;
	W->MeanExcitationEnergy	=	727.0*eV;
	W->ElectronDensity	=	W->CalculateElectronDensity();
	W->PlasmaEnergy		=	W->CalculatePlasmaEnergy();
	db.push_back(W);

	//Lead
	material* Pb = new material();
	Pb->atomic_number	=	82;
	Pb->name		=	"Lead";
	Pb->symbol		=	"Pb";
	Pb->sigma_pN_total	=	2.960;
	Pb->sigma_pN_inelastic	=	1.77;
	Pb->sigma_Rutherford	=	0.00907;
	Pb->dEdx		=	3.40;
	Pb->rho			=	11.35;
	Pb->A			=	207.2;
	Pb->sigma		=	4.8077E6;
	Pb->X0			=	0.005612;
	Pb->density		=	11350;
	Pb->b_N			=	455.3;
	Pb->MeanExcitationEnergy=	823.0*eV;
	Pb->ElectronDensity	=	Pb->CalculateElectronDensity();
	Pb->PlasmaEnergy	=	Pb->CalculatePlasmaEnergy();
	db.push_back(Pb);
}

//Try and find the material we want
material* Material_Database::find_material(string symbol)
{
	//Iterator for moving through the vector of material pointers
	std::vector<material*>::iterator position = db.begin();
	//cout << symbol << endl;

	while(position != db.end())
	{

		if((*position)->symbol == symbol)
		{
			//This should return a pointer to the material we are interested in.
			return *position;
		}
	//Move to the next material	
	*position++;
	}

	//did not find the material, this is bad.
	if(position == db.end())
	{
		if((*position)->symbol == symbol)
		{
			//This should return a pointer to the material we are interested in.
			return *position;
		}
		std::cerr << "Requested aperture material not found. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}
}
