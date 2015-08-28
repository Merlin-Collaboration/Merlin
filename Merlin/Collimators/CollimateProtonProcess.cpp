/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		2010	 RJB
// Modified:	08.01.15 Haroon Rafique		
// Last Edited: 20.07.15 HR
// 
/////////////////////////////////////////////////////////////////////////

#include <iterator>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>

#include "merlin_config.h"

#include "AcceleratorModel/Aperture.h"
#include "AcceleratorModel/Apertures/CollimatorAperture.h"
#include "AcceleratorModel/Apertures/InterpolatedApertures.h"
#include "AcceleratorModel/StdComponent/Collimator.h"

#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

#include "Collimators/CollimateProtonProcess.h"
#include "Collimators/ScatteringProcess.h"
#include "Collimators/ScatteringModel.h"

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"

using namespace std;
using namespace Collimation;

namespace {

	using namespace ParticleTracking;

void OutputIndexParticles(const PSvectorArray lost_p, const list<size_t>& lost_i, ostream& os)
	{
		PSvectorArray::const_iterator p = lost_p.begin();
		list<size_t>::const_iterator ip = lost_i.begin();

		while(p!=lost_p.end())
		{
		    os << std::setw(12) << right << *ip;
		    os << *p;
		    ++p;
		    ++ip;
		}
	}

} // end anonymous namespace

namespace ParticleTracking {

CollimateProtonProcess::CollimateProtonProcess (int priority, int mode, std::ostream* osp)
        : CollimateParticleProcess(priority, mode, osp), dustset(0)
{ 
	
}

bool CollimateProtonProcess::DoScatter(Particle& p)
{  // returns true if particle survives, false if it dies
	//~ std::cout << "\nCollimateProtonProcess::DoScatter" << endl;
	double E0 = currentBunch->GetReferenceMomentum();
	
	Collimator* C = static_cast<Collimator*> (currentComponent); 
	const Aperture *colap = C->GetAperture();
	
	//set scattering model
	SetScatter(scattermodel);
	int smode = scattermodel->GetScatteringPhysicsModel();		

	// this is only ever called for a collimator
//TODO
	//change to step length
	double lengthtogo = GetMaxAllowedStepSize();
	double z = 0;
	int n = 0;

	do{
		double E1 = E0 * (1 + p.dp());
		//Note that pathlength should be calculated with E0
		//~ std::cout << "\nCollimateProtonProcess::DoScatter Call PathLength" << endl;
		double xlen = scattermodel->PathLength(C->p, E0);
		
		double E2 = 0;
		double minLength = min(xlen, lengthtogo);		
		
		double zstep = xlen * sqrt( 1 - p.xp()*p.xp() - p.yp()*p.yp() );
		p.x() += xlen * p.xp();
		p.y() += xlen * p.yp();
		//~ std::cout << "CollimateProtonProcess: p.x() = " << p.x() << std::endl;

//Energy Loss		
		if(smode == 1 || smode == 4){
			//Advanced
			//~ std::cout << "\nCollimateProtonProcess::DoScatter Call Advanced EnergyLoss" << endl;
			scattermodel->EnergyLoss(p, minLength, C->p, E0);
			E2 = E0 * (1 + p.dp());
		}
		else{
			//Simple
			//~ std::cout << "\nCollimateProtonProcess::DoScatter Call Simple EnergyLoss" << endl;
			scattermodel->EnergyLoss(p, minLength, C->p, E0, E1);
			E2 = E0 * (1 + p.dp());
		}
			if(p.dp() < -0.95){
				scattermodel->DeathReport(p, minLength, currentComponent->GetComponentLatticePosition(), lostparticles);
				if(dustset){outputdustbin->Dispose(*currentComponent, xlen, p);}
				return true;}	
//Straggle
		//~ std::cout << "\nCollimateProtonProcess::DoScatter Call Straggle" << endl;
		scattermodel->Straggle(p, minLength, C->p, E1, E2);
			if(colap->PointInside( (p.x()), (p.y()), z+=zstep )){			
			   	//cout << "\n colap->PointInside = " << colap->PointInside( (p.x() + z*p.xp()), (p.y() + z*p.yp()), z ) << " back in aperture"<< endl;
				return false;
			}

//Check if exited Collimator - path length longer than collimator remaining length; not scattered
		if(xlen>lengthtogo) {
			//cout << "Exited Collimator: xlen > LengthToGo" << endl;
			return false;					
		}				
		

//Scattering - use E2 as particle has lost energy up untill this point at which it has scattered
		//~ std::cout << "\nCollimateProtonProcess::DoScatter Call ParticleScatter" << endl;
		if(!scattermodel->ParticleScatter(p, C->p, E2)){			
			//~ p.x() += lengthtogo * p.xp();			//Just for outputting kicks
			//~ p.y() += lengthtogo * p.yp();	
			//~ cout << "CollimateProtonProcess::DoScatter: Inelastic Scatter: particle lost" << endl;	
			scattermodel->DeathReport(p, xlen, currentComponent->GetComponentLatticePosition(), lostparticles);
			//~ cout << "CollimateProtonProcess::DoScatter: Inelastic Scatter: Dispose()" << endl;	
			if(dustset){outputdustbin->Dispose(*currentComponent, xlen, p);}
			return true;
		}
		if(p.dp() < -0.95){
			//cout << "dp <-0.95: particle lost" << endl;
			scattermodel->DeathReport(p, xlen, currentComponent->GetComponentLatticePosition(), lostparticles);
			if(dustset){outputdustbin->Dispose(*currentComponent, xlen, p);}
			return true;
		}
		//cout << "z = " << z << endl;
		//cout << "\nRestart DoScatter n = " << n  << endl;
			
		lengthtogo -= xlen;
		++n;

	} while(true);

	//~ std::cout << "\nCollimateProtonProcess::DoScatter End DoScatter" << endl;
}

/*
//OLD - used when projecting and not tracking protons forward in 10cm bins
bool CollimateProtonProcess::DoScatter(Particle& p)
{  // returns true if particle survives, false if it dies
	double E0 = currentBunch->GetReferenceMomentum();
	
	Collimator* C = static_cast<Collimator*> (currentComponent); 
	const Aperture *colap = C->GetAperture();
	
	//set scattering mode
	CrossSections* SCS = CS_iterator->second;
	int smode = SCS->scat_type;	

	// this is only ever called for a collimator
//TODO
	//change to step length
	double lengthtogo = C->GetLength();
	double z = next_s;
	int n = 0;

	do{
		double E1 = E0 * (1 + p.dp());
		//Note that pathlength should be calculated with E0
		double xlen = C->scatter->PathLength(C->p, E0);
		
		double E2 = 0;
		double minLength = min(xlen, lengthtogo);
		z += minLength;

//Energy Loss		
		if(smode == 1 || smode == 4){
			//Advanced
			E2 = C->scatter->EnergyLoss(p, minLength, C->p, E0, E1);
		}
		else{
		//Simple
		C->scatter->EnergyLoss(p, minLength, C->p, E0, E1);
		E2 = E0 * (1 + p.dp());
		}
			if(p.dp() < -0.95){
				C->scatter->DeathReport(p, minLength, currentComponent->GetComponentLatticePosition(), lostparticles);
				outputdustbin->Dispose(*currentComponent, xlen, p);
				return true;}	
//Straggle
		C->scatter->Straggle(p, minLength, C->p, E1, E2);
			if(colap->PointInside( (p.x() + z*p.xp()), (p.y() + z*p.yp()), z )){			
			   	//cout << "\n colap->PointInside = " << colap->PointInside( (p.x() + z*p.xp()), (p.y() + z*p.yp()), z ) << " back in aperture"<< endl;
				return false;
			}

//Check if exited Collimator - path length longer than collimator remaining length; not scattered
		if(xlen>lengthtogo) {
			//cout << "Exited Collimator: xlen > LengthToGo" << endl;
			return false;					
		}				
		

//Scattering - use E2 as particle has lost energy up untill this point at which it has scattered
		if(!C->scatter->ParticleScatter(p, C->p, E2)){			
			p.x() += lengthtogo * p.xp();			//Just for outputting kicks
			p.y() += lengthtogo * p.yp();	
			//cout << "Inelastic Scatter: particle lost" << endl;	
			C->scatter->DeathReport(p, xlen, currentComponent->GetComponentLatticePosition(), lostparticles);
			outputdustbin->Dispose(*currentComponent, xlen, p);
			return true;
		}
		if(p.dp() < -0.95){
			//cout << "dp <-0.95: particle lost" << endl;
			C->scatter->DeathReport(p, xlen, currentComponent->GetComponentLatticePosition(), lostparticles);
			outputdustbin->Dispose(*currentComponent, xlen, p);
			return true;
		}
		//cout << "z = " << z << endl;
		//cout << "\nRestart DoScatter n = " << n  << endl;
			
		lengthtogo -= xlen;
		++n;

	} while(true);

}*/

void CollimateProtonProcess::SetScatter(ScatteringModel* sm){

//set scattering mode
	if (sm->GetScatteringPhysicsModel() == 0){	//SixTrack
		(sm)->AddProcess(new Inelastic());
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
	}
	else if (sm->GetScatteringPhysicsModel() ==1){ //ST + Adv Ionisation
		(sm)->AddProcess(new Inelastic());
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
	}
	else if (sm->GetScatteringPhysicsModel() == 2){ //ST + Adv Elastic
		(sm)->AddProcess(new Inelastic());
		(sm)->AddProcess(new ElasticpN());
		(sm)->AddProcess(new Elasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());	
		(sm)->AddProcess(new SixTrackRutherford());	
	}
	else if (sm->GetScatteringPhysicsModel() == 3){ //ST + Adv SD
		(sm)->AddProcess(new Inelastic());
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
	}
	else if (sm->GetScatteringPhysicsModel() == 4){ //Merlin
		(sm)->AddProcess(new Inelastic());
		(sm)->AddProcess(new ElasticpN());
		(sm)->AddProcess(new Elasticpn());
		(sm)->AddProcess(new SingleDiffractive());
		(sm)->AddProcess(new Rutherford());
	}
	else{
		std::cout << "Warning ScatteringModel::SetScatterType: No scatter type selected, no ScatteringProcesses set by default - may be set by user" << endl;
	}

}

void CollimateProtonProcess::SetScatteringModel(Collimation::ScatteringModel* s){
	scattermodel = s;
};



}; // end namespace ParticleTracking

