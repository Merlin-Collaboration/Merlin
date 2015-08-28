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
#include "Collimators/DiffractiveScatter.h"
#include "Collimators/ElasticScatter.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

#include "Random/RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace Collimation;

ScatteringModel::ScatteringModel(){}

double ScatteringModel::PathLength(Material* mat, double E0){ 
	
	//~ std::cout << "\nScatteringModel::PathLength: Start PathLength()" << std::endl;
	
	static double lambda; 
	CrossSections* CurrentCS = new CrossSections();
	CrossSections* NewCS;
	
	//~ std::cout << "\nScatteringModel::PathLength: CrossSections object made" << std::endl;
	//~ std::cout << "ScatteringModel::PathLength: map.size = " << stored_cross_sections.size()  << std::endl;
	//~ std::cout << "ScatteringModel::PathLength: mat->GetSymbol() = " << mat->GetSymbol()  << std::endl;
	
	CS_iterator = stored_cross_sections.find(mat->GetSymbol());		
	//~ std::cout << "ScatteringModel::PathLength: CS Iterator sent to find material" << std::endl;
	
	// If find gets to the end of the stored_cross_sections map, there is no value stored
	if (CS_iterator == stored_cross_sections.end() ){
		
		//~ std::cout << "\nScatteringModel::PathLength: Stored cross sections not found " << endl;
	
		//No previously calculated CrossSections, start from scratch	
		NewCS = new CrossSections(mat, E0, ScatteringPhysicsModel);
		
		//~ stored_cross_sections.insert(mat->GetSymbol(), NewCS);		
		stored_cross_sections.insert(std::map< string, Collimation::CrossSections* >::value_type(mat->GetSymbol(), NewCS));
		
		//Set iterator to correct position
		CS_iterator = stored_cross_sections.find(mat->GetSymbol());	
		
		CurrentCS = NewCS;
		if (CurrentCS == NewCS){
			//~ cout << "\n\tScatteringModel::PathLength: CurrentCS == NewCS" << endl;
			//~ delete NewCS;
		}
		else {cout <<  "\n\tWarning: ScatteringModel::PathLength: CurrentCS != NewCS" << endl;}
		
		//~ std::cout << "\nScatteringModel::PathLength: Calculating fractions " << endl;
		//Find fractions of cross sections
		double sigma = 0;
		int i = 0;
		vector<ScatteringProcess*>::iterator p;

		for(p = Processes.begin(); p != Processes.end(); p++){
			(*p)->Configure(mat, CurrentCS);
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
	}
	else{		
		
		//~ std::cout << "\nScatteringModel::PathLength: Stored cross sections found " << endl;
		//Should return a pointer to the CrossSections we require
		CurrentCS = CS_iterator->second;
		
		//Make sure that the CrossSections are for the same case (scattering etc)
		if (CurrentCS == CS_iterator->second){
			//~ cout << "\n\tScatteringModel::PathLength: CurrentCS == NewCS"<< endl;
		}	
		else {
			cout <<  "\n\tWarning: ScatteringModel::PathLength: CurrentCS != StoredCS, recalculating CrossSections" << endl;
			NewCS = new CrossSections(mat, E0, ScatteringPhysicsModel);			
			stored_cross_sections.insert(std::map< string, Collimation::CrossSections* >::value_type(mat->GetSymbol(), NewCS));		
			
			//Set iterator to correct position
			CS_iterator = stored_cross_sections.find(mat->GetSymbol());
			
			CurrentCS = NewCS;
			if (CurrentCS == NewCS){
				//~ cout << "\n\tScatteringModel::PathLength: CurrentCS == NewCS" << endl;
				//~ delete NewCS;
			}
			else {cout <<  "\n\tWarning: ScatteringModel::PathLength: CurrentCS != NewCS" << endl;}			
		}		
		 		
	}
	
	//~ std::cout << "\nScatteringModel::PathLength: calculating lambda " << endl;	
	//Calculate mean free path
	//~ lambda=(mat->GetAtomicNumber())/(sigma*mat->GetDensity()*(1.E-22*Avogadro));
	lambda = CurrentCS->GetTotalMeanFreePath();
	//~ std::cout << "ScatteringModel::PathLength: lambda = " << lambda << endl;
	//~ std::cout << "ScatteringModel::PathLength: returning mean free path " << endl;
	return -(lambda)*log(RandomNG::uniform(0,1));
}             

//Simple energy loss
void ScatteringModel::EnergyLoss(PSvector& p, double x, Material* mat, double E0, double E1){
	double dp = x * (mat->GetSixtrackdEdx());
    p.dp() = ((E1 - dp) - E0) / E0;
}

//Advanced energy loss
void ScatteringModel::EnergyLoss(PSvector& p, double x, Material* mat, double E0){
	
	double E1 = E0 * (1 + p.dp());
	double gamma = E1/(ProtonMassMeV*MeV);
	double beta = sqrt(1 - ( 1 / (gamma*gamma)));
	double I = mat->GetMeanExcitationEnergy()/eV;
	
	double land = RandomNG::landau();
	
	double tmax = (2*ElectronMassMeV * beta * beta * gamma * gamma ) / (1 + (2 * gamma * (ElectronMassMeV/ProtonMassMeV)) + pow((ElectronMassMeV/ProtonMassMeV),2))*MeV;
	
	static const double xi1 = 2.0 * pi * pow(ElectronRadius,2) * ElectronMass * pow(SpeedOfLight,2);
	double xi0 = xi1 * mat->GetElectronDensity();	
	double xi = (xi0 * x /(beta*beta)) / ElectronCharge * (eV/MeV);

	double C = 1 + 2*log(I/mat->GetPlasmaEnergy()/eV);
	double C1 = 0;
	double C0 = 0;

	if((I/eV) < 100)
	{
		if(C <= 3.681)
		{
			C0 = 0.2;
			C1 = 2.0;
		}
		else
		{
			C0 = 0.326*C - 1.0;
			C1 = 2.0;
		}
	}
	else	//I >= 100eV
	{
		if(C <= 5.215)
		{
			C0 = 0.2;
			C1 = 3.0;
		}
		else
		{
			C0 = 0.326*C - 1.5;
			C1 = 3.0;
		}
	}
	double delta = 0;
	
	//Density correction
	double ddx = log10(beta*gamma);
	if(ddx > C1)
	{
		delta = 4.606*ddx - C;
	}
	else if(ddx >= C0 && ddx <= C1)
	{
		double m = 3.0;
		double xa = C /4.606;
		double a = 4.606 * (xa - C0) / pow((C1-C0),m);
		delta = 4.606*ddx -C + a*pow((C1 - ddx),m);
	}
	else
	{
		delta = 0.0;
	}

	double tcut = 2.0*MeV;
	tcut = tmax;
	
	//Mott Correction
	double G = pi*FineStructureConstant*beta/2.0;
	double q = (2*(tmax/MeV)*(ElectronMassMeV) )/(pow((0.843/MeV),2));
	double S = log(1+q);
	double L1 = 0.0;
	double yL2 = FineStructureConstant/beta;

	double L2sum = 1.202001688211;	//Sequence limit calculated with mathematica
	double L2 = -yL2*yL2*L2sum;

	double F = G - S + 2*(L1 + L2);
	double deltaE = xi * (log(2 * ElectronMassMeV * beta*beta * gamma*gamma * (tcut/MeV)/pow(I/MeV,2)) - (beta*beta)*(1 + ((tcut/MeV)/(tmax/MeV))) - delta + F - 1.0 - euler);
	deltaE = xi * (log(2 * ElectronMassMeV * beta*beta * gamma*gamma * xi /pow(I/MeV,2)) - (beta*beta) - delta + F + 0.20);

	double dp = ((xi * land) - deltaE) * MeV;

    p.dp() = ((E1 - dp) - E0) / E0;
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
	//~ std::cout << "ScatteringModel::Straggle: p.x() = " << p.x() << std::endl;
	
 }


bool ScatteringModel::ParticleScatter(PSvector& p, Material* mat, double E0){ 
	//~ std::cout << "\nScatteringModel::ParticleScatter" << std::endl;
	double r = RandomNG::uniform(0,1);
	for(int i = 0; i<fraction.size(); i++)  
	{ 
		//~ std::cout << "ScatteringModel::ParticleScatter: p.x() = " << p.x()<< std::endl;
		//~ std::cout << "\nScatteringModel::ParticleScatter r = " << r << " fraction[i] = " << fraction[i] << std::endl;
	    r -= fraction[i]; 
	    if(r<0)
	    {
	        return Processes[i]->Scatter(p, E0);
	    }
	}
	cout << " should never get this message : \n\tScatteringModel::ParticleScatter : scattering past r < 0" << endl;
	return Processes[0]->Scatter(p, E0);
}

void ScatteringModel::DeathReport(PSvector& p, double x, double position, vector<double>& lost){
    //cout << " particle absorbed\n";
    //~ std::cout << "\nScatteringModel::DeathReport" << std::endl;
	double r = RandomNG::uniform(0,1);
	double pos = x + position;
	lost.push_back(pos);
}
	
void ScatteringModel::SetScatterType(int st){	
	ScatteringPhysicsModel = st;	
}
	
