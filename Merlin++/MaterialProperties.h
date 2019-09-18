/*
 * MaterialProperties.h
 *
 *  Created on: 28 Apr 2017
 *      Author: roger
 */

#ifndef MATERIALPROPERTIES_H_
#define MATERIALPROPERTIES_H_

#include <vector>

/**
 * MaterialProperties.h
 *
 *  Contains the physical quantities for a material
 *  Used for collimation and in wakefields
 *
 *  Created on: 16 Aug 2017
 *      Author: roger
 */

#include <iostream>
#include <string>
#include <map>

using namespace std;

class MaterialProperties
{
public:
	double Z, A, density, dEdx, sigma_R, sigma_I, sigma_T, X0;
	double dEdxGeVpermetre, lambda;
	map<string, double> *extra;

	virtual double A_R()
	{
		return A;
	}
	virtual double A_H()
	{
		return A;
	}
	MaterialProperties()
	{
		lambda = dEdxGeVpermetre = A = density = dEdx = sigma_R = sigma_I = sigma_T = Z = X0 = 0;
		extra = nullptr;
	}
	// keep child class happy
	MaterialProperties(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
	void Update();
	MaterialProperties* EnergyScale(double E);
	virtual ~MaterialProperties();
	MaterialProperties(const MaterialProperties&); // copy constructor
	MaterialProperties& operator=(const MaterialProperties&); // copy assignment
	void SetExtra(string, ...);
	double GetExtra(string);
	bool HaveExtra(string);
};

class Mixture: public MaterialProperties
{
public:
	int N;
	vector<double> V_R;
	vector<double> V_H;
	vector<double> V_A;
	double A_H(); // random nucleus for elastic scattering
	double A_R(); // random nucleus for Rutherford scattering
	Mixture(map<string, MaterialProperties*> dict, vector<string> names, vector<double>
		proportions, double d);
};

ostream& operator<<(ostream& o, MaterialProperties& d);

#endif /* MATERIALPROPERTIES_H_ */
