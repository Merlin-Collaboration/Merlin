#include "Collimators/Material_Database.hpp"
#include "Collimators/Material.hpp"

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

Material_Database::Material_Database() : db(0) //size of the db vector, aka number of materials, increase as required.
{
	//Here we create new materials, add their properties, then push them onto a vector for manipulation.
	//Need to get references for all data.
/*
References:

atomic_number	
name		
symbol		
sigma_pN_total
sigma_pN_inelastic
sigma_Rutherford	
dEdx	
rho	
A	
sigma
X0
density		
Particle data group: http://pdg.lbl.gov/2009/AtomicNuclearProperties/
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
	//Be->A			=	9.012182;
	Be->A			=	9.01;		//SIXTRACK EDIT
	Be->sigma		=	3.08E7;
	//Be->X0		=	0.3528;
	Be->X0			=	0.353;		//SIXTRACK EDIT
	Be->density		=	1848;
	Be->b_N			=	74.7;
	db.push_back(Be);

	//Carbon (graphite)
	material* C = new material();
	C->atomic_number	=	6;
	C->name			=	"Carbon";
	C->symbol		=	"C";
	C->sigma_pN_total	=	0.331;
	C->sigma_pN_inelastic	=	0.231;
	C->sigma_Rutherford	=	0.000076;
//	C->dEdx			=	0.68;
	C->dEdx			=	0.75;		//SIXTRACK EDIT
//	C->rho			=	2.265;		//Check
	C->rho			=	2.26;		//SIXTRACK EDIT		//4.52 for C2?
//	C->A			=	12.0107;
	C->A			=	12.01;		//SIXTRACK EDIT
	C->sigma		=	7.14E4;
	//C->X0			=	18.8E-2;	//Not for graphite
	C->X0			=	0.1932;
	C->density		=	2210;
	C->b_N			=	70.0;
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
//	Al->rho			=	1.204;		//wrong? - yes!
	Al->rho			=	2.70;		//SIXTRACK EDIT
//	Al->A			=	26.981538;
	Al->A			=	26.98;		//SIXTRACK EDIT
	Al->sigma		=	35.64E6;
//	Al->X0			=	0.08897;
	Al->X0			=	0.089;		//SIXTRACK EDIT
	Al->density		=	2699;
	Al->b_N			=	120.3;
	db.push_back(Al);

	//Copper
	material* Cu = new material();
	Cu->atomic_number	=	29;
	Cu->name		=	"Copper";
	Cu->symbol		=	"Cu";
	Cu->sigma_pN_total	=	1.232;
	Cu->sigma_pN_inelastic	=	0.782;
	Cu->sigma_Rutherford	=	0.00153;
	Cu->dEdx		=	2.69;
	Cu->rho			=	8.96;
//	Cu->A			=	63.546;
	Cu->A			=	63.55;		//SIXTRACK EDIT
	Cu->sigma		=	5.98E7;
//	Cu->X0			=	0.01436;
	Cu->X0			=	0.0143;		//SIXTRACK EDIT
	Cu->density		=	8960;
	Cu->b_N			=	217.8;
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
//	W->A			=	183.84;
	W->A			=	183.85;		//SIXTRACK EDIT
	W->sigma		=	0.177E4;
//	W->X0			=	0.003504;
	W->X0			=	0.0035;		//SIXTRACK EDIT
	W->density		=	19300;
	W->b_N			=	440.3;
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
//	Pb->A			=	207.2;
	Pb->A			=	207.19;		//SIXTRACK EDIT
	Pb->sigma		=	4.8077E6;
//	Pb->X0			=	0.005612;
	Pb->X0			=	0.0056;		//SIXTRACK EDIT
	Pb->density		=	11350;
	Pb->b_N			=	455.3;
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
