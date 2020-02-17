#include "MaterialData.h"
#include "MaterialProperties.h"
#include "PhysicalUnits.h"
#include <iomanip>
#include <vector>
#include <stdarg.h>

using namespace PhysicalUnits;

MaterialData::~MaterialData()  // destructor
{
	map<string, MaterialProperties*>::iterator p = property.begin();
	while(p != property.end())
	{
		delete p->second; p++;
	}
}

StandardMaterialData::StandardMaterialData()   // constructor
{
	string names[] = {"Al", "Be", "C", "Cu", "Pb", "Si", "W", "O", "Fe", "Ni", "Mo"};
	double ST[] = {0.634, 0.268, 0.331, 1.232, 2.96, 0.653, 2.767, -1, -1, -1, -1};
	double SI[] = {0.421, 0.199, 0.231, 0.782, 1.71, 0.432, 1.65, -1, -1, -1, -1};
	double A[] = {26.9815, 9.01218, 12.011, 63.546, 207.2, 28.1, 183.84, 15.9994, 55.845, 58.6934, 95.96};
	double SR[] = {0.00034, 0.000035, 0.000076, 0.00153, 0.00907, 0.00039, 0.00768, 0.0133766, -1, -1, -1};
	double dE[] = {1.724, 1.594, 1.745, 1.403, 1.123, 1.664, 1.145, -1, -1, -1, -1};
	double X[] = {24.01, 65.19, 42.7, 12.86, 6.37, 21.82, 6.76, -1, -1, -1, -1};
	double rho[] = {1.204, 1.848, 2.265, 8.96, 11.35, 2.33, 19.3, 1.14, 7.879, 8.9, 10.2};  // g/cc
	double Z[] = {13., 4., 6., 29., 82., 14., 74., 8, 26, 28, 42};
	int n = sizeof(names) / sizeof(string);
	for(int i = 0; i < n; i++)
	{
		property[names[i]] = new MaterialProperties(A[i], ST[i], SI[i], SR[i], dE[i], X[i], rho[i] * gram / cc, Z[i]);
	}

	string names2[] = {"Al", "Be", "C", "Cu", "Pb", "Si", "W", "Fe", "Ni", "Mo"};
	double conductivity[] = {3.564E7, 3.08E7, 71400, 5.98E7, 4.8077E6, 1, 1770, 1.04E7, 1.44E7, 1.87E7};
	for(int i = 0; i < int(sizeof(names2) / sizeof(string)); i++)
		property[names2[i]]->SetExtra("conductivity", conductivity[i]);

	string names4[] = {"Be", "C", "Al", "Cu", "W", "Pb"};
	double sixtracknucslope[] = {74.7, 70, 120.3, 217.8, 440.3, 455.3};
	for(int i = 0; i < int(sizeof(names4) / sizeof(string)); i++)
		property[names4[i]]->SetExtra("sixtracknuclearslope", sixtracknucslope[i]);

	string names5[] = {"Be", "C", "O", "Al", "Fe", "Ni", "Cu", "Mo", "W", "Pb"};
	double MeanExcitationEnergy[] = {63.7, 78.0, 95.0, 166.0, 286.0, 311.0, 322.0, 424.0, 727.0, 823.0 }; //  in eV
	for(int i = 0; i < int(sizeof(names5) / sizeof(string)); i++)
		property[names5[i]]->SetExtra("MeanExcitationEnergy", MeanExcitationEnergy[i] * eV);

	property["AC150K"] = new MaterialProperties(*property["C"]);
	property["AC150K"]->dEdx = 0.68;
	property["AC150K"]->density = 1.650 * gram / cc;

	MakeMixture("IT180", "W Ni Cu", 0.95, 0.35, 0.15, 18.06 * gram / cc);

	MakeMixtureByWeight("Glidcop", "Cu Al O", 0.9972, .0015, .0013, 8.93 * gram / cc);
	property["Glidcop"]->SetExtra("conductivity", 5.38E7);

}

void StandardMaterialData::UseSixTrackValues()
{
	string names3[] = {"Be", "C", "Al", "Fe", "Cu", "W", "Pb"};
	double sixtracktotnu[] = {0.268, .331, 1.232, 2.767, 2.960};
	for(int i = 0; i < int(sizeof(names3) / sizeof(string)); i++)
		property[names3[i]]->sigma_T = sixtracktotnu[i];

	string names6[] = {"C", "Al", "Cu", "W"};
	double sixtrackinelasticnuclear[] = {0.231, 0.421, 0.782, 1.65};
	for(int i = 0; i < int(sizeof(names6) / sizeof(string)); i++)
		property[names6[i]]->sigma_I = sixtrackinelasticnuclear[i];

	string names7[] = {"C", "O", "Al", "Cu", "W"};
	double sixtrackrutherford[] = {0.000076, 0.0133766, 0.00034, 0.00153, 0.00768};
	for(int i = 0; i < int(sizeof(names7) / sizeof(string)); i++)
		property[names7[i]]->sigma_R = sixtrackrutherford[i];

	string names8[] = {"C", "Al", "Cu", "W", "AC150K"};
	double sixtrackdedx[] = {0.68, 0.81, 2.69, 5.79, 0.68};
	for(int i = 0; i < int(sizeof(names8) / sizeof(string)); i++)
		property[names8[i]]->dEdx = sixtrackdedx[i];
}

void MaterialData::PrintTable()
{
	cout << "  Name    A    Z    sig T  sig I  sig R   dEdx  density Rad length\n";
	for(map<string, MaterialProperties*>::iterator it = property.begin(); it != property.end(); it++)
		cout << setw(7) << it->first << " " << *(it->second) << endl;
}
void MaterialData::MakeMixture(string name, string s, ...)
{
	vector<string> elements;
	vector<double> proportions;
	va_list numbers;
	va_start(numbers, s);
	unsigned long int i;
	do
	{
		i = s.find(" ");
		string want = s.substr(0, i);
		elements.push_back(want);
		proportions.push_back(va_arg(numbers, double));
		s = s.substr(i + 1);
	} while(i != string::npos);
	double density = va_arg(numbers, double);
	// cout << " density is " << density << endl;
	va_end(numbers);
	property[name] = new Mixture(property, elements, proportions, density);
}
void MaterialData::MakeMixtureByWeight(string name, string s, ...)
{
	vector<string> elements;
	vector<double> proportions;
	va_list numbers;
	va_start(numbers, s);
	unsigned long int i;
	do
	{
		i = s.find(" ");
		string want = s.substr(0, i);
		elements.push_back(want);
		proportions.push_back(va_arg(numbers, double) / property[want]->A);
		s = s.substr(i + 1);
	} while(i != string::npos);
	double density = va_arg(numbers, double);
	va_end(numbers);
	property[name] = new Mixture(property, elements, proportions, density);
}

ostream& operator<<(ostream& s, MaterialData* M)
{
	s << "  Name    A    Z    sig T  sig I  sig R   dEdx  density Rad length\n";
	for(map<string, MaterialProperties*>::iterator it = (M->property).begin(); it != M->property.end(); it++)
		cout << setw(7) << it->first << " " << *(it->second) << endl;

	return s;
}
