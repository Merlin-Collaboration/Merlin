/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		08.01.15 Haroon Rafique	
// Modified:	08.01.15 Haroon Rafique		
// Last Edited: 20.07.15 HR
// 
/////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>

#include "AcceleratorModel/AcceleratorComponent.h"

#include "Collimators/Dustbin.h"

#include "Exception/MerlinException.h"

using namespace std;
namespace ParticleTracking {

Dustbin::Dustbin(OutputType ot)
{
	otype = ot;
}

LossMapDustbin::LossMapDustbin(OutputType ot)
{
	otype = ot;
}

void LossMapDustbin::Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle)
{

	//~ cout << "\nLossMapDustbin called" << endl;

	if (currentComponent != &currcomponent)
	{	
		currentComponent = &currcomponent;
	}
	temp.reset();

	temp.ElementName = currentComponent->GetQualifiedName().c_str();
	temp.s = currentComponent->GetComponentLatticePosition();
	// pos is the lost position within the element, position is the exact loss position in the lattice
	temp.position = (pos + temp.s);
	temp.length = currentComponent->GetLength();
	temp.lost = 1;
		
	//calculate 10cm interval - move to 10cm binning?
	double inter = 0.0;
	bool fin = 0;

	do
	{
		if ( (pos >= inter) && (pos < (inter+0.1)) )
		{
			temp.interval = inter;
			fin = 1;
		}
		else 
		{
			inter += 0.1;
		}
	}
	while(fin == 0);
	
	//obsolete
	//if(currentComponent->GetType() =="Collimator"){temp.temperature = 0;}				
	//else if(currentComponent->GetType() =="SectorBend"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="Quadrupole"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="SkewQuadrupole"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="XYKicker"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="XCor"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="YCor"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="BPM"){temp.temperature = 1;}
	//else if(currentComponent->GetType() =="Drift"){temp.temperature = 1;}
	//else {temp.temperature = 2;}	

	//For LHC Loss Maps
	if(currentComponent->GetType() =="Collimator"){temp.temperature = 0;}				
	else if(temp.ElementName.substr(0, 14) =="SectorBend.MBW"){temp.temperature = 2;}	
	else if(temp.ElementName.substr(0, 15) =="SectorBend.MBXW"){temp.temperature = 2;}
	else if(temp.ElementName.substr(0, 14) =="Quadrupole.MQW"){temp.temperature = 2;}
	else if(temp.ElementName.substr(0, 9) =="XCor.MCBW"){temp.temperature = 2;}
	else if(temp.ElementName.substr(0, 9) =="YCor.MCBW"){temp.temperature = 2;}
	else if(temp.ElementName.substr(0, 9) =="BPM.BPMWE"){temp.temperature = 2;}
	else if(temp.ElementName.substr(0, 10) =="Drift.MCBW"){temp.temperature = 2;}
	else if( (temp.s > 1.9790884399999788E+04) && (temp.s <= 2.0244198399999801E+04 ) && (temp.ElementName.substr(0, 11) == "Drift.DRIFT") ) {temp.temperature = 2;} //Drifts in IR7
	else {temp.temperature = 1;}	

	temp.p = particle;
	//pushback vector
	DeadParticles.push_back(temp);

}

void LossMapDustbin::Finalise()
{	
	//First sort DeadParticles according to s
	sort(DeadParticles.begin(), DeadParticles.end(), Compare_LossData); 	

	cout << "DUSTBIN:: DeadParticles.size() = " <<  DeadParticles.size() << endl;

	int outit = 0;
	int total = 0;

	switch(otype) 
	{
	case nearestelement:
		for(vector<LossData>::iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			// Start at s = min and push back the first LossData
			if (OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it); 
			}				
			// If old element ++loss
			if (it->ElementName == OutputLosses[outit].ElementName)
			{
				OutputLosses[outit].lost +=1;
			}
			// If new element OutputLosses.push_back
			else
			{
				OutputLosses.push_back(*it);	
				outit++;				
			}
		}
	break;
	case precise:
		for(vector<LossData>::iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			// Start at s = min and push back the first LossData
			if (OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it); 
			}	
			// If position is equal
			if (it->position == OutputLosses[outit].position)
			{
				OutputLosses[outit].lost +=1;
			}
			// If new element outit.push_back
			else
			{
				OutputLosses.push_back(*it);	
				outit++;				
			}
		}
	break;
	case tencm:
		for(vector<LossData>::iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			//if no losses are yet stored
			if (OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it); 
			}	
			// If in the same bin ++loss
			if ((it->ElementName == OutputLosses[outit].ElementName) && (it->interval == OutputLosses[outit].interval))
			{
				OutputLosses[outit].lost +=1;
			}
			// If new element outit.push_back and set loss to 1
			else
			{
				OutputLosses.push_back(*it);	
				outit++;				
				OutputLosses[outit].lost =1;
			}
		}
	break;
	};
	cout << "DUSTBIN:: OutputLosses.size() = " << OutputLosses.size() << endl;
	cout << "DUSTBIN:: Total losses = " << total << endl;
}

void LossMapDustbin::Output(std::ostream* os)
{
	switch(otype) 
	{	
		case nearestelement:
		for(vector <LossData>::iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << setw(34) << left << (*its).ElementName;
			(*os) << setw(34) << left << (*its).s;
			(*os) << setw(16) << left << (*its).lost;
			(*os) << setw(16) << left << (*its).temperature;
			(*os) << setw(16) << left << (*its).length;
			(*os) << endl;
		}
		break;

		case precise:
		for(vector <LossData>::iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << setw(34) << left << (*its).ElementName;
			(*os) << setw(34) << left << (*its).s;
			(*os) << setw(34) << left << (*its).position;
			(*os) << setw(16) << left << (*its).lost;
			(*os) << setw(16) << left << (*its).temperature;
			(*os) << endl;

		}
		break;

		case tencm:
		for(vector <LossData>::iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << setw(34) << left << (*its).ElementName;
			(*os) << setw(34) << left << (*its).s;
			(*os) << setw(16) << left << (*its).interval;
			(*os) << setw(16) << left << (*its).lost;
			(*os) << setw(16) << left << (*its).temperature;
			(*os) << endl;

		}
		break;
	}

}

}; // End namespace ParticleTracking
