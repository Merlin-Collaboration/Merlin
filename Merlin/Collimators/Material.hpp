#ifndef _MATERIAL_HPP_
#define _MATERIAL_HPP_

#include <string>

using namespace std;

class material //(particle.type)?
{

public:

	size_t atomic_number;		//Atomic number
	string name;			//Element/compound name
	string symbol;			//Elemental symbol: what input file will check for
	double sig_pN_tot;		//proton nucleus total cross section =sigma_el_pn + sigma_el_pN + sigma_SD_pn + sigma_inel_pN
	double sig_pN_in;		//proton nucleus inelastic cross section (sigma_inel_pN)
	double sig_R;			//Rutherford scattering cross section (sigma_R)
	double dEdx;			//
	double rho;			//Material Density, g/cm^3
	double A;			//Atomic mass
	double sigma;			//electrical conductivity (sigma) = 1/electrical resisitivity (Ohm*m)e-1
	double X0;			//Radiation Length??
	double density;			//Material density, kg/m^3

//Sub-structure for scattering?
//foo.electron_scatter.sigma
//foo.proton_scatter.sigma
//psedocode:
/*
if(particle == "proton")
{
	scatter = &proton_scatter;
}
else if(particle == "electron")
{
	scatter = &electron_scatter;
}
else
{
	std::cout << "No appropriate scattering data found" << std::endl;
	exit(EXIT_FAILURE);
}
*/

//But if we have electron/hadron collisions, both will be required
//But each collimator will only see one type of beam

private:

};

#endif
