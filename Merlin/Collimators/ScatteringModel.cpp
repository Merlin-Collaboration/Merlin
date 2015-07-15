/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		
// Modified:	
// Last Edited: 18.03.15 HR
// 
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>

#include "Collimators/ScatteringModel.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalConstants.h"

#include "Random/RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

ScatteringModel::ScatteringModel(){lastmat=0;}

double ScatteringModel::PathLength(Material* mat, double E0){ 
	static double lambda; 
	if(mat != lastmat){ // recalculate cross sections if material has changed
		lastmat = mat;
		cout << " recomputing cross sections \n";
		double sigma = 0;
		int i = 0;
		vector<ScatteringProcess*>::iterator p;

		for(p = Processes.begin(); p != Processes.end(); p++){
			(*p)->Configure(mat, E0);
			fraction[i]= (*p)->sigma;
			cout << (*p)->GetProcessType() << "\t\t sigma = " << (*p)->sigma << " barns" << endl;
			sigma+= fraction[i];
			++i;
		}

		for(int j=0;j<fraction.size();j++){
			cout << " Process " << j << " total sigma " << setw(10) << setprecision(4) << sigma << "barns";
			fraction[j] /= sigma;
			cout << " fraction " << setw(10) << setprecision(4) << fraction[j] << endl;
		}  

		lambda=(mat->GetAtomicNumber())/(sigma*mat->GetDensity()*(1.E-22*Avogadro));
		cout << " Lambda = " << scientific << lambda << endl;
	}
	return -(lambda)*log(RandomNG::uniform(0,1));
	//return -(2*lambda)*log(RandomNG::uniform(0,1));
}             


double ScatteringModel::EnergyLoss(PSvector& p, double x, Material* mat, double E0, double E1){
     p.dp() -=  (x * (mat->GetRadiationLengthInM()) )/E1; 
     return (E0 * (1 + p.dp()));
}


//HR 29Aug13
void ScatteringModel::Straggle(PSvector& p, double x, Material* mat, double E0, double E2){
		
	// working
	static const double root12 = sqrt(12.0);
	double scaledx=x/mat->GetRadiationLengthInM();
	double Eav = (E0+E2)/2.0;
	//sixtrack plus
	double theta0 = 13.6*MeV * sqrt (scaledx) * (1.0 + 0.038 * log(scaledx)) / Eav; 
	//sixtrack
	//double theta0 = 13.6*MeV * sqrt (scaledx) / Eav; 
	
	//JM
	double theta_plane_x = RandomNG::normal(0,1) * theta0;
	double theta_plane_y = RandomNG::normal(0,1) * theta0;
	
	double x_plane = RandomNG::normal(0,1) * x * (theta0/root12) + x * theta_plane_x/2;
	double y_plane = RandomNG::normal(0,1) * x * (theta0/root12) + x * theta_plane_y/2;
	
	p.x ()  += x_plane;
	p.xp () += theta_plane_x; 
	p.y ()  += y_plane;
	p.yp () += theta_plane_y; 
	
 }


bool ScatteringModel::ParticleScatter(PSvector& p, Material* mat, double E0){ 
	double r = RandomNG::uniform(0,1);
	for(int i = 0; i<fraction.size(); i++)  
	{ 
	    r -= fraction[i]; 
	    if(r<0)
	    {
	        return Processes[i]->scatter(p, E0);
	    }
	}
	cout << " should never get this message : \n\tScatteringModel::ParticleScatter : scattering past r < 0" << endl;
	return Processes[0]->scatter(p, E0);
}

void ScatteringModel::DeathReport(PSvector& p, double x, double position, vector<double>& lost){
    //cout << " particle absorbed\n";
	double pos = x + position;
	lost.push_back(pos);
	
}




