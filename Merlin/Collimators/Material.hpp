#ifndef _MATERIAL_HPP_
#define _MATERIAL_HPP_

#include <string>

using namespace std;

class material
{
public:
//selected Sixtrack variable names in quotes

	size_t atomic_number;		//Atomic number
	string name;			//Element/compound name
	string symbol;			//Elemental symbol: what input file will check for
	double sigma_pN_total;		//proton nucleus total cross section =sigma_el_pn + sigma_el_pN + sigma_SD_pn + sigma_inel_pN
	double sigma_pN_inelastic;	//proton nucleus inelastic cross section (sigma_inel_pN)
	double sigma_Rutherford;	//Rutherford scattering cross section (sigma_R)
	double dEdx;			// "dpodx"
	double rho;			//Material Density, g/cm^3
	double A;			//Atomic mass
	double sigma;			//electrical conductivity (sigma) = 1/electrical resisitivity (Ohm*m)e-1
	double X0;			//Radiation Length
	double density;			//Material density, kg/m^3
	double b_N;			//Nuclear slope
	double ElectronDensity;		//Electron density: calculated from other input data
	double ElectronCriticalEnergy;
	double MeanExcitationEnergy;
	double PlasmaEnergy;
	
	double CalculateElectronDensity();
	double CalculatePlasmaEnergy();
	double CalculateMeanExcitationEnergy();
};

#endif
