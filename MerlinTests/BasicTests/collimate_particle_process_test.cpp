#include "../tests.h"
#include <iostream>


#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/Components.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchTypes.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "Collimators/CollimateParticleProcess.h"

/* Create a bunch of particle, and check that the correct number survive various sized apertures
 *
 */


using namespace std;
using namespace PhysicalUnits;

int main(int argc, char* argv[])
{
	
	AcceleratorModelConstructor* ctor = new AcceleratorModelConstructor();
	ctor->NewModel();

	AcceleratorComponent *drift = new Drift("d1",1*meter);

	RectangularAperture *rect_app = new RectangularAperture(21*millimeter, 100*millimeter);
	ctor->AppendComponent(*drift);

	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;

	double beam_energy = 7000.0;
	ProtonBunch* myBunch = new ProtonBunch(beam_energy,1);

	double coords[][4] = {
		{0,0,0,0},
		{10,0,10,0},
		{-20,0,20,0},
		{10,0,0,0},
		{20,0,0,0},
		{0,0,10,0},
		{0,0,20,0},
		{0,0,-30,0},
	};
	size_t npart = sizeof(coords) / sizeof(coords[0]);
	cout << "npart " << npart << endl;
	
	vector<Particle> pcoords;
	for(size_t i=0; i<npart; i++){
		Particle p(0);
		p.x() = coords[i][0]*millimeter;
		p.xp() = coords[i][1];
		p.y() = coords[i][2]*millimeter;
		p.yp() = coords[i][3];
		pcoords.push_back(p);
		//myBunch->AddParticle(p);
	}

	AcceleratorModel::RingIterator ring = theModel->GetRing();
	ParticleTracker* tracker = new ParticleTracker(ring,myBunch);
	
	        //#define COLL_AT_ENTRANCE 1
	        //#define COLL_AT_CENTER 2
	        //#define COLL_AT_EXIT 4
	CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(2,4);
	tracker->AddProcess(myCollimateProcess);
	

	
	drift->SetAperture(rect_app);


	ProtonBunch* test_bunch;
	vector<Particle> pcoords2;
	//
	// Test a few rectangular appertures with tracking
	//
	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy,1, pcoords2);
	rect_app->SetFullWidth(70*millimeter);
	rect_app->SetFullHeight(70*millimeter);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 8);
	delete test_bunch;


	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy,1, pcoords2);
	rect_app->SetFullWidth(70*millimeter);
	rect_app->SetFullHeight(50*millimeter);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 7);
	delete test_bunch;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy,1, pcoords2);
	rect_app->SetFullWidth(50*millimeter);
	rect_app->SetFullHeight(70*millimeter);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 8);
	delete test_bunch;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy,1, pcoords2);
	rect_app->SetFullWidth(30*millimeter);
	rect_app->SetFullHeight(30*millimeter);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 4);
	delete test_bunch;

	delete rect_app;

	delete myBunch;
	delete tracker;
	delete theModel;

}

