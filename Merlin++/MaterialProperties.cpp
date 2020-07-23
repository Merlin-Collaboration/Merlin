/*
 * MaterialProperties.cpp
 *
 *  Created on: 16 Aug 2017
 *      Author: roger
 */

#include "MaterialProperties.h"

#include <iomanip>
#include <stdarg.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include "RandomNG.h"
#include "PhysicalUnits.h"
#include "MerlinException.h"

using namespace std;
using namespace PhysicalUnits;

MaterialProperties::MaterialProperties(double p1, double p2, double p3, double p4, double p5, double p6, double p7,
	double p8)
{
	A = p1;
	sigma_T = p2;
	sigma_I = p3;
	sigma_R = p4;
	dEdx = p5;
	X0 = p6;
	density = p7;
	Z = p8;
	extra = new map<string, double>; // will be deleted by destructor
	Update();
}

void MaterialProperties::Update()
{
	lambda = A * 1.E-3 / ((sigma_T + sigma_R) * density * 1.E-28 * 6.022E23);
}

MaterialProperties* MaterialProperties::EnergyScale(double E)
{
	MaterialProperties* p = new MaterialProperties(A, sigma_T, sigma_I, sigma_R, dEdx, X0, density, Z);
	// adjust dEdx and cross sections with new energy
	p->Update();
	return p;
}

MaterialProperties::MaterialProperties(const MaterialProperties& a)
{ //copy constructor
	density = a.density;
	A = a.A;
	sigma_T = a.sigma_T;
	sigma_I = a.sigma_I;
	sigma_R = a.sigma_R;
	dEdx = a.dEdx;
	X0 = a.X0;
	Z = a.Z;
	extra = new  map<string, double>(*a.extra); // deep copy. Deleted by dtor
}

MaterialProperties& MaterialProperties::operator=(const MaterialProperties& a)
{ // copy assignment
	if(this != &a)
	{ // beware of self-assignment
		density = a.density;
		A = a.A;
		sigma_T = a.sigma_T;
		sigma_I = a.sigma_I;
		sigma_R = a.sigma_R;
		dEdx = a.dEdx;
		X0 = a.X0;
		Z = a.Z;
		extra = new  map<string, double>(*a.extra);   //deep copy. deleted by dtor
	}
	return *this;
}

ostream& operator<<(ostream& o, MaterialProperties& d)
{
	o << setiosflags(ios::fixed) << setw(5) << setprecision(1) << d.A
	  << setw(5) << setprecision(1) << d.Z
	  << setw(7) << setprecision(3) << d.sigma_T
	  << setw(7) << d.sigma_I << setw(9) << setprecision(5) << d.sigma_R << setw(6) << setprecision(3) << d.dEdx
	  << setw(7) << d.density / (gram / cc)
	  << " " << fixed
	  << setw(7) << setprecision(1) << d.X0;
	if(d.extra && d.extra->size() > 0)
	{
		for(map<string, double>::iterator i = d.extra->begin(); i != d.
			extra->end(); i++)
			o << " " << i->first << " " << std::scientific << setw(9) << setprecision(2) << i->second << std::fixed;
	}
	return o;
}
//
void MaterialProperties::SetExtra(string props, ...)
{
	/**
	 * Avoids the user having to mess with brackets and contents-of
	 * and can handle an arbitrary number of arguments
	 * though they must match the space-delimited list in the string
	 *  e.g.  SetExtra("one two",1.0,2.0);
	 */
	va_list arguments;
	va_start(arguments, props);
	istringstream s(props);
	string prop;
	while(s >> prop)
	{
		(*extra)[prop] = va_arg(arguments, double);
	}
	va_end(arguments);
}

double MaterialProperties::GetExtra(string prop)
{
	return (*extra)[prop];
}

bool MaterialProperties::HaveExtra(string prop)
{
	map<string, double>::iterator it;
	it = extra->find(prop);
	return it != extra->end();
}

MaterialProperties::~MaterialProperties()
{
	if(extra)
		delete extra;
	extra = nullptr;
}

double Mixture::A_H()  // random nucleus for elastic scattering
{
	float r = RandomNG::uniform(0., 1.);
	for(int i = 0; i < N - 1; i++)
	{
		if((r -= V_H[i]) <= 0)
			return V_A[i];
	}
	return V_A[N - 1];
}

double Mixture::A_R()  // random nucleus for Rutherford scattering
{
	float r = RandomNG::uniform(0., 1.);
	for(int i = 0; i < N - 1; i++)
	{
		if((r -= V_R[i]) <= 0)
			return V_A[i];
	}
	return V_A[N - 1];
}

Mixture::Mixture(map<string, MaterialProperties*> dict, vector<string> names, vector<double> proportions, double d)
{
	N = proportions.size();
	if(names.size() != proportions.size())
	{
		throw MerlinException("Mixture::Mixture() length of names and proportions do not match");
	}
	float suma = 0, sum = 0, sumH = 0, sumR = 0;
	for(int i = 0; i < N; i++)
	{
		if(dict.count(names[i]) == 0)
		{
			throw MerlinException("Mixture::Mixture() component not found: " + names[i]);
		}
		sum += proportions[i];
		suma += proportions[i] * dict[names[i]]->A;
		V_A.push_back(dict[names[i]]->A);
		V_R.push_back(dict[names[i]]->sigma_R);
		sumR += dict[names[i]]->sigma_R;
		V_H.push_back(dict[names[i]]->sigma_T); // MUST CHECK. WANT NUCLEAR ELASTIC
		sumH += dict[names[i]]->sigma_T;
	}
	Z = A = sigma_T = sigma_I = sigma_R = dEdx = 0;
	float RX0 = 0;
	for(int i = 0; i < N; i++)
	{
		Z += dict[names[i]]->Z * proportions[i] / sum;
		A += dict[names[i]]->A * proportions[i] / sum;
		sigma_T += dict[names[i]]->sigma_T * proportions[i] / sum;
		sigma_R += dict[names[i]]->sigma_R * proportions[i] / sum;
		sigma_I += dict[names[i]]->sigma_I * proportions[i] / sum;
//    dE/dx and (reciprocal) Radiation length use mass fractions rather than num erical fractions
		//     as they will be multiplied by the density before use
		dEdx += dict[names[i]]->dEdx * proportions[i] * dict[names[i]]->A / suma
		;
		RX0 += (1. / dict[names[i]]->X0) * proportions[i] * dict[names[i]]->A / suma;
		V_R[i] /= sumR;
		V_H[i] /= sumH;
	}
	X0 = 1 / RX0;
	density = d;
	// cout << "  Mixture mixture densioty " << d << endl;
	extra = new map<string, double>; // will be deleted by destructor
	Update();
}

/*
   ostream& operator<<(ostream& o, MaterialProperties d)
   {
    o << setiosflags(ios::fixed)
      << setw(5) << setprecision(1) << d.A
      << setw(5) << setprecision(1) << d.Z
      << setw(7) << setprecision(3) << d.sigma_T << setw(7) << d.sigma_I
      << setw(9) << setprecision(5) << d.sigma_R
      << setw(6) << setprecision(3) << d.dEdx << setw(7) << d.density
      << " "
      << fixed << setw(7) << setprecision(1) << d.X0;
    return o;
   }
 */
