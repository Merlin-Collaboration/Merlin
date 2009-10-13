// ILC Spin Rotator Example
// -----------------------------------------------------------------------
// Original implemention: Andy Wolski
// modified by          : Nick Walker 5/10/2005
//
// This application performs particle tracking through a possible ILC
// spin rotator. Spin tracking is also included.
//
// Routines are also provided for studying tolerances (sensitivities) of
// magnets (in this case quadrupoles)
//
// The example also shows the use of a user-defined SimulationOuput class
// for controlling the recording of results.
//

#include "merlin_config.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "Random/RandomNG.h"
#include "MADInterface/MADInterface.h"
#include "BeamModel/BeamData.h"
#include "AcceleratorModel/StdComponent/StandardMultipoles.h"
#include "AcceleratorModel/Frames/TComponentFrame.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/SpinParticleProcess.h"
#include "BeamDynamics/ParticleTracking/SynchRadParticleProcess.h"
#include "BasicTransport/NormalTransform.h"

#define BEAMENERGY 5.0*GeV

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

// Forward global function declarations
// (definitions can be found after main()

SpinParticleBunch* ConstructSpinParticleBunch(const BeamData&, const SpinVector&, int nparts);

// Output class declaration
// Here we define a class which we will use to output information that we need

class SpinTrackingOutput : public SimulationOutput {
public:
	SpinTrackingOutput(const std::string& fname);
	std::ostream& os() { return fos; }
protected:
	void Record(const ComponentFrame* frame, const Bunch* bunch);
	void RecordInitialBunch(const Bunch* bunch) { Output("INITIAL",*bunch);}
	void RecordFinalBunch(const Bunch* bunch)   { Output("FINAL",*bunch); }

private:
	void Output(const std::string& label, const Bunch& bunch);
	std::ofstream fos;
};

// ---------------------------------------------------------------------
// MAIN programme
// ---------------------------------------------------------------------

int main()
{
	// Initialise random number generator
	RandomNG::init(123456);

	// Simulation flags
	const bool inc_synch_rad = false;

	try {
		// Construct the AcceleratorModel
		// from a lattice file produced by MAD
		MADInterface madi("../lattices/SpinRotator.lattice.txt", BEAMENERGY);

		ofstream madlog("mad.log");
		madi.SetLogFile(madlog);
		madi.SetLoggingOn();
		// If we are modelling SR, scale the magnets for
		// the mean energy loss
		madi.ScaleForSynchRad(inc_synch_rad);

		AcceleratorModel* theModel = madi.ConstructModel();

		// Define the initial beam parameters
		BeamData beam;

		double gamma = (BEAMENERGY/MeV)/ElectronMassMeV;
		beam.p0      = BEAMENERGY;
		beam.emit_x  = 8.00e-6/gamma;
		beam.emit_y  = 0.02e-6/gamma;
		beam.sig_z   = 6.00e-3;
		beam.sig_dp  = 1.50e-3;

		beam.beta_x  = 20.012342132996;
		beam.beta_y  = 20.012435560926;

		// Construct a bunch of 10^4 particles, where each particle has vertical spin.
		SpinParticleBunch* spinBunch0 = ConstructSpinParticleBunch(beam,SpinVector(0,1,0),10000);

		// Construct a ParticleTracker to perform the tracking
		// including the spin tracking process
		ParticleTracker tracker(theModel->GetBeamline());
		tracker.AddProcess(new SpinParticleProcess(1));

		if(inc_synch_rad){
			SynchRadParticleProcess* srp = new SynchRadParticleProcess(2,true);
			srp->IncludeQuadRadiation(false); // only model SR in dipoles
			srp->SetNumComponentSteps(1);     // make one step though elements
			tracker.AddProcess(srp);
		}

		// Set the output module for the simulation
		SpinTrackingOutput* spOut = new SpinTrackingOutput("spin-results.dat");
		spOut->output_initial = true;
		spOut->output_final = true;
		spOut->output_all = true;
		// spOut->AddIdentifier("Quadrupole.*"); // Output at exit of quadrupoles
		tracker.SetOutput(spOut);

		// Do the tracking. Note we track a copy of our original bunch
		SpinParticleBunch* spinBunch = new SpinParticleBunch(*spinBunch0);
		tracker.Track(spinBunch);

		// Write the phase space co-ordinates
		// of each particle in the bunch to a file
		// output format is (1 parlicle/line)
		// ct0 E0 x x' y y' dct dE/E px py pz 
//		ofstream trackingLog("SpinTracking.dat");
//		spinBunch->Output(trackingLog);

		delete spinBunch;
		tracker.SetOutput(0);
		delete spOut;

		// Quadrupole sensitivities
		// Having set up the tracker and the output, we will now
		// see how the vertical motion of the quadrupoles affect
		// the results.
		
		// Again we use a SpinTrackingOuput object, but this
		// time we only want to have the final result per quadrupole
		// moved.
		spOut = new SpinTrackingOutput("quad-alignment-sensitivity.dat");
		spOut->output_initial = false;
		spOut->output_final = true;
		spOut->output_all = false;
		tracker.SetOutput(spOut);

		std::vector< TComponentFrame<Quadrupole>* > quads;
		theModel->ExtractTypedComponents(quads);
		cout<<quads.size()<<" Quadrupoles in lattice"<<endl;

		const double dy = 1.0*micrometer;

		for(size_t i = 0; i<quads.size(); i++) {
			string quadID = (quads[i]->GetComponent()).GetQualifiedName();
			cout<<"  adjusting "<<quadID<<"..."<<flush;
			
			(spOut->os())<<setw(20)<<left<<quadID;
			
			quads[i]->TranslateY(dy);
			spinBunch = new SpinParticleBunch(*spinBunch0);
			tracker.Track(spinBunch);
			delete spinBunch;
			quads[i]->ClearLocalFrameTransform();

			cout<<"done"<<endl;
		}

		// Clean up and quit
		delete theModel;
		cout<<"Finished!"<<endl;
	
	} catch(MerlinException& error) {
		cerr<<error.Msg()<<endl;
		return -99;
	}
	return 0;
}

// Global function definitions

SpinParticleBunch* ConstructSpinParticleBunch(const BeamData& beam, const SpinVector& spin, int nparts)
{
	ParticleBunchConstructor pbc(beam,nparts,normalDistribution);
	ParticleBunch* aBunch = pbc.ConstructParticleBunch();
	SpinParticleBunch* spinBunch = new SpinParticleBunch(BEAMENERGY);

	for(PSvectorArray::iterator it = aBunch->begin(); it!=aBunch->end(); it++) {
		PSvector p=*it;
		spinBunch->AddParticle(p,spin);
	}

	delete aBunch; // No longer needed.
	return spinBunch;
}

// class SpinTrackingOutput definitions

SpinTrackingOutput::SpinTrackingOutput(const std::string& fname)
:fos(fname.c_str())
{
	if(!fos)
		throw MerlinException(string("Probem openning output file: ")+fname);
}

void SpinTrackingOutput::Record(const ComponentFrame* frame, const Bunch* bunch)
{
	Output((*frame).GetComponent().GetQualifiedName(),*bunch);
}

// simple macro to output floating point number
#define WRITE_OS(os,w,p,data) (os)<<scientific<<setw(w)<<setprecision(p)<<(data)

void SpinTrackingOutput::Output(const std::string& label, const Bunch& bunch)
{
	// static_cast here should be OK as we are only
	// dealing with SpinParticleBunch (should use dynamic_cast to be completely
	// safe).
	const SpinParticleBunch& spinBunch = static_cast<const SpinParticleBunch&>(bunch);

	PSmoments S;
	spinBunch.GetMoments(S);

	double p0 = spinBunch.GetReferenceMomentum();
	p0*=1+S.mean(ps_DP);

	double gamma = p0/MeV/ElectronMassMeV;
	double gex = gamma*ProjectedEmittance(S,ps_X,ps_XP);
	double gey = gamma*ProjectedEmittance(S,ps_Y,ps_YP);
	double z = spinBunch.GetReferenceTime();

	// Remove the <ydp> and <y'dp> correlations and 
	// re-calculate the emittance if the energy
	// spread is not zero.
	double geyc = gey;
	if(S.var(ps_DP) != 0) {
		double s36 = S(ps_Y,ps_DP);
		double s46 = S(ps_YP,ps_DP);
		double dp2 = S.var(ps_DP);
		double s33 = S.var(ps_Y)-s36*s36/dp2;
		double s34 = S(ps_Y,ps_YP)-s36*s46/dp2;
		double s44 = S.var(ps_YP)-s46*s46/dp2;
		geyc = gamma*sqrt(s33*s44-s34*s34);
	}

	SpinVector pa = spinBunch.GetAverageSpin();
	double spinMag = sqrt(pa.x()*pa.x()+pa.y()*pa.y()+pa.z()*pa.z());

	fos<<setw(20)<<left<<label<<right;	// element label
	WRITE_OS(fos,12,3,z);				// beamline location (m)
	WRITE_OS(fos,12,3,p0);				// mean momentum (GeV/c)
	WRITE_OS(fos,12,3,gey);				// normalised vertical emittance (m)
	WRITE_OS(fos,12,3,geyc);			//      "        "         "      "  dispersion-corrected 
	WRITE_OS(fos,12,3,gex);				// normalised horizontal emittance (m)
	WRITE_OS(fos,12,3,S.mean(ps_Y));    // vertical beam centroid displacement (m)
	WRITE_OS(fos,12,3,S.std(ps_Y));     // vertical RMS beam size (m)
	WRITE_OS(fos,12,3,pa.x());          // average spin vector components
	WRITE_OS(fos,12,3,pa.y());
	WRITE_OS(fos,12,3,pa.z());
	WRITE_OS(fos,12,3,spinMag);         // average spin vector magnetude
	fos<<endl;
}

