#include "NumericalUtils/PhysicalUnits.h"
#include "MADInterface/MADInterface.h"
#include "AcceleratorModel/Supports/MagnetMover.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"

#include "BPMVectorBuffer.h"

#define BEAMENERGY 5.0*GeV

typedef vector< MagnetMover* > MagnetMoverList;

using namespace PhysicalUnits;
using namespace ParticleTracking;

int main()
{
	// Construct the AcceleratorModel from a lattice file produced by MAD.
	MADInterface madi("../lattices/MERLINFodo.lattice.txt", BEAMENERGY);

	ofstream madlog("mad.log");
	madi.SetLogFile(madlog);
	madi.SetLoggingOn();

	AcceleratorModel* theModel = madi.ConstructModel();


	// Extract a list of magnet movers from the AcceleratorModel
	// and translate the 20th mover 20 microns vertically.
	MagnetMoverList magnetMovers;
	theModel->ExtractTypedElements(magnetMovers);
	magnetMovers[20]->SetY(20.0e-6);


	// Construct a bunch of particles to track through the lattice.
	// Here we just add a single particle on the reference trajectory.
	ParticleBunch* theBunch = new ParticleBunch(BEAMENERGY);

	PSvector p(0);
	theBunch->AddParticle(p);


	// Construct a ParticleTracker to perform the tracking.
	ParticleTracker tracker(theModel->GetBeamline(), theBunch);


	// Construct a BPMBuffer to record the bunch centroid at each BPM.
	BPMVectorBuffer* bpmVecBuffer = new BPMVectorBuffer();
	BPM::SetDefaultBuffer(bpmVecBuffer);

	// Do the tracking.
	tracker.Run();


	// Write the tracking results to a file.
	ofstream bpmLog("Trajectory.dat");
	vector<BPM::Data>& theBPMBuffer = bpmVecBuffer->BPMReading;
	for(vector<BPM::Data>::iterator bpm_iter=theBPMBuffer.begin(); bpm_iter!=theBPMBuffer.end(); bpm_iter++) {
		bpmLog<<std::setw(14)<<(bpm_iter->x).value;
		bpmLog<<std::setw(14)<<(bpm_iter->x).error;
		bpmLog<<std::setw(14)<<(bpm_iter->y).value;
		bpmLog<<std::setw(14)<<(bpm_iter->y).error;
		bpmLog<<endl;
	};


	BPM::SetDefaultBuffer(0);
	delete bpmVecBuffer;
	
	delete theBunch;
	delete theModel;

	cout<<"Finished!"<<endl;
	return 0;
}
