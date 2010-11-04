#ifndef _MATERIAL_HPP_
#define _MATERIAL_HPP_

#include <string>

using namespace std;

class material
{
public:
//selected Sixtrack variable names in quotes

	size_t atomic_number;		//Atomic number	"zatom"
	string name;			//Element/compound name
	string symbol;			//Elemental symbol: what input file will check for
	double sigma_pN_total;		//proton nucleus total cross section =sigma_el_pn + sigma_el_pN + sigma_SD_pn + sigma_inel_pN
	double sigma_pN_inelastic;	//proton nucleus inelastic cross section (sigma_inel_pN)
	double sigma_Rutherford;	//Rutherford scattering cross section (sigma_R)
	double dEdx;			// "dpodx"
	double rho;			//Material Density, g/cm^3 "rho"
	double A;			//Atomic mass "anuc"
	double sigma;			//electrical conductivity (sigma) = 1/electrical resisitivity (Ohm*m)e-1
	double X0;			//Radiation Length "radl"
	double density;			//Material density, kg/m^3
	double b_N;			//Nuclear slope
};

#endif
