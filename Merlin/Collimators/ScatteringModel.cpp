/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		
// Modified:	07.09.15 Haroon Rafique		
// Last Edited: 07.09.15 HR
// 
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstring>

#include "Collimators/ScatteringModel.h"
#include "Collimators/DiffractiveScatter.h"
#include "Collimators/ElasticScatter.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

#include "Random/RandomNG.h"

using namespace ParticleTracking;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace Collimation;
using namespace std;

ScatteringModel::ScatteringModel()
{
	//~ ScatteringModelTurn = 1;	
	ScatterPlot_on = 0;
	JawImpact_on = 0;
}

double ScatteringModel::PathLength(Material* mat, double E0){ 
	
	//~ std::cout << "\nScatteringModel::PathLength: Start PathLength()" << std::endl;
	
	static double lambda; 
	CrossSections* CurrentCS;
	
	//~ std::cout << "\nScatteringModel::PathLength: CrossSections object made" << std::endl;
	//~ std::cout << "ScatteringModel::PathLength: map.size = " << stored_cross_sections.size()  << std::endl;
	//~ std::cout << "ScatteringModel::PathLength: mat->GetSymbol() = " << mat->GetSymbol()  << std::endl;
	
	CS_iterator = stored_cross_sections.find(mat->GetSymbol());		
	//~ std::cout << "ScatteringModel::PathLength: CS Iterator sent to find material" << std::endl;
	
	// If find gets to the end of the stored_cross_sections map, there is no value stored
	if (CS_iterator == stored_cross_sections.end() ){
		
		//~ std::cout << "\nScatteringModel::PathLength: Stored cross sections not found " << endl;
	
		//No previously calculated CrossSections, start from scratch	
		CurrentCS = new CrossSections(mat, E0, ScatteringPhysicsModel);
		
		//~ stored_cross_sections.insert(mat->GetSymbol(), NewCS);		
		stored_cross_sections.insert(std::map< string, Collimation::CrossSections* >::value_type(mat->GetSymbol(), CurrentCS));
		
		//Set iterator to correct position
		CS_iterator = stored_cross_sections.find(mat->GetSymbol());			
		
		//~ std::cout << "\nScatteringModel::PathLength: Calculating fractions " << endl;
		//Find fractions of cross sections
		double sigma = 0;
		int i = 0;
		vector<ScatteringProcess*>::iterator p;

		std::cout << "ScatteringModel::PathLength: MATERIAL = " << mat->GetSymbol()  << std::endl;
		for(p = Processes.begin(); p != Processes.end(); p++){
			(*p)->Configure(mat, CurrentCS);
			fraction[i] = (*p)->sigma;
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
			CurrentCS = new CrossSections(mat, E0, ScatteringPhysicsModel);			
			stored_cross_sections.insert(std::map< string, Collimation::CrossSections* >::value_type(mat->GetSymbol(), CurrentCS));		
			
			//Set iterator to correct position
			CS_iterator = stored_cross_sections.find(mat->GetSymbol());				
		}		
		 		
	}
	
	//Calculate mean free path
	lambda = CurrentCS->GetTotalMeanFreePath();
	//~ std::cout << "ScatteringModel::PathLength: lambda = " << lambda << endl;
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

	double C = 1 + 2*log(I/(mat->GetPlasmaEnergy()/eV));
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
void ScatteringModel::Straggle(PSvector& p, double x, Material* mat, double E1, double E2){

	static const double root12 = sqrt(12.0);
	double scaledx=x/mat->GetRadiationLengthInM();
	double Eav = (E1+E2)/2.0;
	double theta0 = 13.6*MeV * sqrt (scaledx) * (1.0 + 0.038 * log(scaledx)) / Eav; 
	
	double theta_plane_x = RandomNG::normal(0,1) * theta0;
	double theta_plane_y = RandomNG::normal(0,1) * theta0;
	
	double x_plane = RandomNG::normal(0,1) * x * theta0 / root12 + x * theta_plane_x / 2;
	double y_plane = RandomNG::normal(0,1) * x * theta0 / root12 + x * theta_plane_y / 2;
	
	p.x ()  += x_plane;
	p.xp () += theta_plane_x; 
	p.y ()  += y_plane;
	p.yp () += theta_plane_y; 
	
 }


bool ScatteringModel::ParticleScatter(PSvector& p, Material* mat, double E){ 

	double r = RandomNG::uniform(0,1);
	for(int i = 0; i<fraction.size(); i++)  
	{ 
	    r -= fraction[i]; 
	    if(r<0)
	    {
	        return Processes[i]->Scatter(p, E);
	    }
	}
	
	cout << " should never get this message : \n\tScatteringModel::ParticleScatter : scattering past r < 0, r = " << r << endl;

	exit(EXIT_FAILURE);
}

void ScatteringModel::DeathReport(PSvector& p, double x, double position, vector<double>& lost){
	double pos = x + position;
	lost.push_back(pos);
}
	
void ScatteringModel::SetScatterType(int st){	
	ScatteringPhysicsModel = st;	
}

void ScatteringModel::ScatterPlot(Particle& p, double z, int turn, string name){
	
	ScatterPlotData* temp = new ScatterPlotData;
	(*temp).ID = p.id();
	(*temp).x = p.x();
	(*temp).xp = p.xp();
	(*temp).y = p.y();
	(*temp).yp = p.yp();
	(*temp).z = z;
	(*temp).turn = turn;	
	(*temp).name = name;
	
	StoredScatterPlotData.push_back(temp);
}

void ScatteringModel::JawImpact(Particle& p, int turn, string name){
	
	JawImpactData* temp = new JawImpactData;
	(*temp).ID = p.id();
	(*temp).x = p.x();
	(*temp).xp = p.xp();
	(*temp).y = p.y();
	(*temp).yp = p.yp();
	(*temp).ct = p.ct();
	(*temp).dp = p.dp();
	(*temp).turn = turn;	
	(*temp).name = name;
	
	StoredJawImpactData.push_back(temp);
}
	
	
void ScatteringModel::SetScatterPlot(string name, int single_turn){
	ScatterPlotNames.push_back(name);
	ScatterPlot_on = 1;
}

void ScatteringModel::SetJawImpact(string name, int single_turn){
	JawImpactNames.push_back(name);
	JawImpact_on = 1;
}

//~ void ScatteringModel::OutputScatterPlot(std::ostream* os){
void ScatteringModel::OutputScatterPlot(string directory, int seed){
	
	for(vector<string>::iterator name = ScatterPlotNames.begin(); name != ScatterPlotNames.end(); ++name){
		
		std::ostringstream scatter_plot_file;
		scatter_plot_file << directory << "scatter_plot_" << (*name) << "_" << seed << ".txt";
		//~ std::ostringstream scatter_plot_file = directory + "Scatter.txt";
		std::ofstream* os = new std::ofstream(scatter_plot_file.str().c_str());	
		if(!os->good())    {
			std::cerr << "ScatteringModel::OutputJawImpact: Could not open ScatterPlot file for collimator " << (*name) << std::endl;
			exit(EXIT_FAILURE);
		} 
		
		//~ (*os) << "#\tparticle_id\tx\tx'\ty\ty'\tct\tdpctturn" << endl;
		(*os) << "#\tparticle_id\tz\ty\tturn" << endl;
	
		for(vector <ScatterPlotData*>::iterator its = StoredScatterPlotData.begin(); its != StoredScatterPlotData.end(); ++its)
		{
			if( (*its)->name == (*name) ){
				(*os) << setw(10) << setprecision(10) << left << (*its)->ID;
				(*os) << setw(30) << setprecision(20) << left << (*its)->z;
				(*os) << setw(30)<< setprecision(20) << left << (*its)->x;
				(*os) << setw(30) << setprecision(20) << left << (*its)->y;
				//~ (*os) << setw(10)<< setprecision(20) << left << (*its)->ID;
				//~ (*os) << setw(20)<< setprecision(20) << left << (*its)->xp;
				//~ (*os) << setw(20)<< setprecision(20) << left << (*its)->y;
				//~ (*os) << setw(20)<< setprecision(20) << left << (*its)->yp;
				//~ (*os) << setw(20)<< setprecision(20) << left << (*its)->z;
				(*os) << setw(10)<< setprecision(10) << left << (*its)->turn;
				(*os) << endl;
			}	
		}
	}
	
	StoredScatterPlotData.clear();
}

//~ void ScatteringModel::OutputJawImpact(std::ostream* os){
void ScatteringModel::OutputJawImpact(string directory, int seed){
	
	for(vector<string>::iterator name = JawImpactNames.begin(); name != JawImpactNames.end(); ++name){
		std::ostringstream jaw_impact_file;
		jaw_impact_file << directory << "jaw_impact_" << (*name) << "_" << seed << ".txt";
		ofstream* os = new ofstream(jaw_impact_file.str().c_str());	
		if(!os->good())    {
			std::cerr << "ScatteringModel::OutputJawImpact: Could not open JawImpact file for collimator " << (*name) << std::endl;
			exit(EXIT_FAILURE);
		} 
		
		(*os) << "#\tparticle_id\tx\tx'\ty\ty'\tct\tdpctturn" << endl;
	
		for(vector <JawImpactData*>::iterator its = StoredJawImpactData.begin(); its != StoredJawImpactData.end(); ++its)
		{
			if( (*its)->name == (*name) ){
				(*os) << setw(10) << left << setprecision(10) <<  (*its)->ID;
				(*os) << setw(30) << left << setprecision(20) << (*its)->x;
				(*os) << setw(30) << left << setprecision(20) <<  (*its)->xp;
				(*os) << setw(30) << left << setprecision(20) <<  (*its)->y;
				(*os) << setw(30) << left << setprecision(20) <<  (*its)->yp;
				(*os) << setw(30) << left << setprecision(20) << (*its)->ct;
				(*os) << setw(30) << left << setprecision(20) << (*its)->dp;
				(*os) << setw(10) << left << setprecision(10) <<  (*its)->turn;
				(*os) << endl;
			}
		}
	}	
		
	StoredJawImpactData.clear();	
}
