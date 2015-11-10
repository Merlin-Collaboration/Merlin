/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 5.01 (2015)
// 
// Copyright: see Merlin/copyright.txt
//
// Created:		2010	 RJB
// Modified:	07.09.15 Haroon Rafique		
// Last Edited: 03.11.15 HR
// 
/////////////////////////////////////////////////////////////////////////

#include <iterator>
#include <string>
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
        //~ : CollimateParticleProcess(priority, mode, osp), dustset(0), scattermodel(NULL)
        : CollimateParticleProcess(priority, mode, osp), scattermodel(NULL)
{ 
	
}

bool CollimateProtonProcess::DoScatter(Particle& p)
{  // returns true if particle survives, false if it dies
	
	double P0 = currentBunch->GetReferenceMomentum();	
	double E0 = sqrt(P0*P0 + pow(PhysicalConstants::ProtonMassMeV*PhysicalUnits::MeV,2));
	
	//~ cout << "ColProPro: Turn = " << ColParProTurn << endl;
	
	bool scatter_plot = 0;
	bool jaw_impact = 0;
	
	// Length of the collimator
	double coll_length = currentComponent->GetLength();

	double z = currentBunch->int_s;
	double lengthtogo = s-z;		
	
	//~ cout << "\n\t\tColPROPro z = " << z << endl;
	//~ cout << "\t\tColPROPro s = " << s << endl;
	//~ cout << "\t\tColPROPro lengthtogo = " << lengthtogo << endl;
	
	Collimator* C = static_cast<Collimator*> (currentComponent); 
	
	string ColName = currentComponent->GetName();
	
	if( (scattermodel->ScatterPlot_on) && (scattermodel->ScatterPlotName == ColName) ){
		scatter_plot = 1;
		//~ cout << "ColProPro: ScatterPlot ON: ColName = " << ColName << ", ScatterPlotName = " << scattermodel->ScatterPlotName << endl;
	}
	if( (scattermodel->JawImpact_on) && (scattermodel->JawImpactName == ColName) ){
		jaw_impact = 1;
		//~ cout << "ColProPro: JawImpact ON: ColName = " << ColName << ", JawImpactName = " << scattermodel->JawImpactName << endl;
	}
	
	const Aperture *colap = C->GetAperture();
		
	//set scattering model
	if (scattermodel == NULL){
		std::cout << "\nCollimateProtonProcess::SoScatter::WARNING: no ScatteringModel set." << std::endl;
		std::cout << "Use 'myCollimateProcess->SetScatteringModel(myScatter);'" << std::endl;
		std::cout << "Don't forget to assign a scattering mode using 'myScatter->SetScatterType(n);'" << std::endl;
		std::cout << "Where n: 0=SixTrack, 1=ST+Adv. Ionisation, 2=ST+Adv. Elastic, 3=ST+Adv. Single Diffractive, 4=Merlin" << std::endl;
		exit(EXIT_FAILURE);
	}

	int smode = scattermodel->GetScatteringPhysicsModel();		

	while(lengthtogo>0){
		double E1 = E0 * (1 + p.dp());
		//Note that pathlength should be calculated with E0
		
		double xlen = scattermodel->PathLength(C->p, E0);
				
		double E2 = 0;
				
		bool interacted = ( lengthtogo > xlen );
		double step_size = interacted ? xlen : lengthtogo;
		
		double zstep = step_size * sqrt( 1 - p.xp()*p.xp() - p.yp()*p.yp() );
		
		//~ cout << "\t\tColPROPro step_size = " << step_size << endl;
		//~ cout << "\t\tColPROPro z+step_size = " << z+step_size << endl;
		
		
		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
		
//Jaw Impact
		if(jaw_impact && z == 0){
				scattermodel->JawImpact(p);
		}
		
//Scatter Plot
		if(scatter_plot && z == 0){
		//~ if(scatter_plot){
				scattermodel->ScatterPlot(p, z);
		}
		
//Energy Loss		
		if(smode == 1 || smode == 4){
			//Advanced
			scattermodel->EnergyLoss(p, step_size, C->p, E0);
		}
		else{
			//Simple
			scattermodel->EnergyLoss(p, step_size, C->p, E0, E1);			
		}
		
		E2 = E0 * (1 + p.dp());
		
		//~ if(p.dp() < ((1/E0) - 1)){
		if(E2 <=1.0){
			p.ct() = z;
			scattermodel->DeathReport(p, step_size, currentComponent->GetComponentLatticePosition(), lostparticles);
			//~ if(dustset){outputdustbin->Dispose(*currentComponent, (lengthtogo - step_size), p, ColParProTurn);}
			if(dustset){				
				for(DustbinIterator = DustbinVector.begin(); DustbinIterator != DustbinVector.end(); ++DustbinIterator){					
						(*DustbinIterator)->Dispose(*currentComponent, (z+zstep), p, ColParProTurn);
				}			
			}
			return true;
		}	
//MCS
		scattermodel->Straggle(p, step_size, C->p, E1, E2);
		
		if( (E2 < (E0 / 100.0)) ) {
			return false;					
		}	


//Check if (returned to aperture) OR (travelled through length) 
		z+=zstep;
		if(scatter_plot){scattermodel->ScatterPlot(p, z);}
		
		if( (colap->PointInside( (p.x()), (p.y()), z)) || (xlen>lengthtogo) ) {
			//disabling agrees with assmann test case for 0.1m bins
			//~ p.x() += p.xp()*lengthtogo;
			//~ p.y() += p.yp()*lengthtogo;			
			return false;					
		}						

//Scattering - use E2
		if(interacted){
			if(!scattermodel->ParticleScatter(p, C->p, E2)){		
				p.ct() = z;
				scattermodel->DeathReport(p, step_size, currentComponent->GetComponentLatticePosition(), lostparticles);
				//~ if(dustset){outputdustbin->Dispose(*currentComponent, (lengthtogo - step_size), p, ColParProTurn);}
				if(dustset){					
					for(DustbinIterator = DustbinVector.begin(); DustbinIterator != DustbinVector.end(); ++DustbinIterator){					
						(*DustbinIterator)->Dispose(*currentComponent, (z+zstep), p, ColParProTurn);
					}					
				}
				return true;
			}
		}
		
		if( (p.dp() < -0.95) || (p.dp() < -1) ){
			p.ct() = z;
			scattermodel->DeathReport(p, step_size, currentComponent->GetComponentLatticePosition(), lostparticles);
			//~ if(scatter_plot){scattermodel->ScatterPlot(p, z);}
			//~ if(dustset){outputdustbin->Dispose(*currentComponent, (lengthtogo - step_size), p, ColParProTurn);}if(dustset){					
			if(dustset){					
				for(DustbinIterator = DustbinVector.begin(); DustbinIterator != DustbinVector.end(); ++DustbinIterator){					
					(*DustbinIterator)->Dispose(*currentComponent, (z+zstep), p, ColParProTurn);
				}					
			}
			return true;
		}
		
		lengthtogo -= step_size;

//Scatter Plot
		//~ if(scatter_plot){
				//~ scattermodel->ScatterPlot(p, z);
		//~ }		
	}
}

void CollimateProtonProcess::SetScatter(ScatteringModel* sm){
// set scattering mode
// Note that inelastic should always be the last added process
	if (sm->GetScatteringPhysicsModel() == 0){	//SixTrack
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
		(sm)->AddProcess(new Inelastic());
	}
	else if (sm->GetScatteringPhysicsModel() ==1){ //ST + Adv Ionisation	
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
		(sm)->AddProcess(new Inelastic());
	}
	else if (sm->GetScatteringPhysicsModel() == 2){ //ST + Adv Elastic
		(sm)->AddProcess(new ElasticpN());
		(sm)->AddProcess(new Elasticpn());
		(sm)->AddProcess(new SixTrackSingleDiffractive());	
		(sm)->AddProcess(new SixTrackRutherford());			
		(sm)->AddProcess(new Inelastic());
	}
	else if (sm->GetScatteringPhysicsModel() == 3){ //ST + Adv SD
		(sm)->AddProcess(new SixTrackElasticpN());
		(sm)->AddProcess(new SixTrackElasticpn());
		(sm)->AddProcess(new SingleDiffractive());
		(sm)->AddProcess(new SixTrackRutherford());
		(sm)->AddProcess(new Inelastic());
	}
	else if (sm->GetScatteringPhysicsModel() == 4){ //Merlin
		(sm)->AddProcess(new ElasticpN());
		(sm)->AddProcess(new Elasticpn());
		(sm)->AddProcess(new SingleDiffractive());
		(sm)->AddProcess(new Rutherford());
		(sm)->AddProcess(new Inelastic());
	}
	else{
		std::cout << "Warning ScatteringModel::SetScatterType: No scatter type selected, no ScatteringProcesses set by default - may be set by user" << endl;
	}

}

void CollimateProtonProcess::SetScatteringModel(Collimation::ScatteringModel* s){
	scattermodel = s;
	SetScatter(scattermodel);
};

}; // end namespace ParticleTracking

