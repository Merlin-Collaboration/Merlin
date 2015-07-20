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

#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"

using namespace std;

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
        : CollimateParticleProcess(priority, mode, osp) 
{ 
	
}
/*
bool CollimateProtonProcess::DoScatter(Particle& p)
{  // returns true if particle survives, false if it dies
	double E0 = currentBunch->GetReferenceMomentum();
	
	Collimator* C = static_cast<Collimator*> (currentComponent); 
	const Aperture *colap = C->GetAperture();

	// this is only ever called for a collimator
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
		E2 = C->scatter->EnergyLoss(p, minLength, C->p, E0, E1);
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

}
*/
}; // end namespace ParticleTracking

