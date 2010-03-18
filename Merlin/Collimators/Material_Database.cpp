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


	//Beryllium
	material* Be = new material();
	Be->atomic_number	=	4;
	Be->name		=	"Beryllium";
	Be->symbol		=	"Be";
	Be->sig_pN_tot		=	0.268;
	Be->sig_pN_in		=	0.199;
	Be->sig_R		=	0.000035;
	//Be->dEdx		=	1.594;
	Be->dEdx		=	0.55;
	Be->rho			=	1.848;
	Be->A			=	9.012182;
	Be->sigma		=	3.08E7;
	Be->X0			=	35.28E-2;
	db.push_back(Be);

	//Carbon
	material* C = new material();
	C->atomic_number	=	6;
	C->name			=	"Carbon";
	C->symbol		=	"C";
	C->sig_pN_tot		=	0.331;
	C->sig_pN_in		=	0.231;
	C->sig_R		=	0.000076;
	//C->dEdx			=	1.745;
	C->dEdx			=	0.68;
	C->rho			=	2.265;
	C->A			=	12.0107;
	C->sigma		=	7.14E4;
	C->X0			=	18.8E-2;
	db.push_back(C);

	//Aluminium
	material* Al = new material();
	Al->atomic_number	=	13;
	Al->name		=	"Aluminium";
	Al->symbol		=	"Al";
	Al->sig_pN_tot		=	0.634;
	Al->sig_pN_in		=	0.421;
	Al->sig_R		=	0.00034;
	//Al->dEdx		=	1.724;
	Al->dEdx		=	0.81;
	Al->rho			=	1.204;
	Al->A			=	26.981538;
	Al->sigma		=	35.64E6;
	Al->X0			=	8.90E-2;
	db.push_back(Al);

	//Copper
	material* Cu = new material();
	Cu->atomic_number	=	29;
	Cu->name		=	"Copper";
	Cu->symbol		=	"Cu";
	Cu->sig_pN_tot		=	1.232;
	Cu->sig_pN_in		=	0.782;
	Cu->sig_R		=	0.00153;
	//Cu->dEdx		=	1.430;
	Cu->dEdx		=	2.69;
	Cu->rho			=	8.96;
	Cu->A			=	63.546;
	Cu->sigma		=	5.98E7;
	Cu->X0			=	1.44E-2;
	db.push_back(Cu);

	//Tungsten
	material* W = new material();
	W->atomic_number	=	74;
	W->name			=	"Tungsten";
	W->symbol		=	"W";
	W->sig_pN_tot		=	2.767;
	W->sig_pN_in		=	1.65;
	W->sig_R		=	0.00768;
	//W->dEdx			=	1.145;
	W->dEdx			=	5.79;
	W->rho			=	19.3;
	W->A			=	183.84;
	W->sigma		=	0.177E4;
	W->X0			=	0.351E-2;
	db.push_back(W);

	//Lead
	material* Pb = new material();
	Pb->atomic_number	=	82;
	Pb->name		=	"Lead";
	Pb->symbol		=	"Pb";
	Pb->sig_pN_tot		=	2.960;
	Pb->sig_pN_in		=	1.77;
	Pb->sig_R		=	0.00907;
	//Pb->dEdx		=	1.123;
	Pb->dEdx		=	3.40;
	Pb->rho			=	11.35;
	Pb->A			=	207.2;
	Pb->sigma		=	4.8077E6;
	Pb->X0			=	0.562E-2;
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
		std::cerr << "Requested aperture material not found. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}
}
