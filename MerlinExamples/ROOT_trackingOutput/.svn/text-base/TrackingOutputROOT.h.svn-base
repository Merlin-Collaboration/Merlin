#ifndef _h_TrackingOutputROOT
#define _h_TrackingOutputROOT

#include "BeamDynamics/TrackingSimulation.h"
#include "utility/StringPattern.h"
#include "BeamModel/BeamData.h"

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/SMPTracking/SMPBunch.h"

#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"

using namespace ParticleTracking;
using namespace SMPTracking;

//--------------------------------------------------------------------
//
// Usage examples:
//
// A tracking tree
// 1. one tree, automatic ROOT file creation (merlin.root) 
//    TrackingOutputROOT* to = new TrackingOutputROOT("mytree")
//
// 2. many trees, automatic file name (merlin.root)
//    TrackingOutputROOT* to = new TrackingOutputROOT()
//                    to->NewTree("name1")
//                           ...
//                    to->NewTree("name2")
//
// 3. one tree, explicitly given ROOT file 
//    TFile* rfile = new TFile("mymerlin.root")
//    TrackingOutputROOT* to = new TrackingOutputROOT("mytree",rfile)
//
// 4. many trees, explicitly given ROOT file 
//    TFile* rfile = new TFile("mymerlin.root")
//    TrackingOutputROOT* to = new TrackingOutputROOT("",rfile)
//                    to->NewTree("name1")
//                           ...
//                    to->NewTree("name2")
//
//............................................................................................
//
// B bunch tree
// 1.   At a certain position 
//	to->DumpBunchAt("INITIAL");         // unique position by id
//	to->DumpBunchAt("Quadrupole.*");    // at any(!) quad - pattern
//      to->DoBunch(false);                 // switch bunch tree (temporarily) off
//      to->DoBunch();                      // switch on - default is on
//
// 2. A specific, given bunch
//	ParticleBunch* bunch = ParticleBunchConstructor(*beam,1000).ConstructParticleBunch();
//	to->DumpBunch(bunch,"bunch_before_tracking");
                
//--------------------------------------------------------------------


class TrackingOutputROOT : public SimulationOutput {
public:
	TrackingOutputROOT(const std::string& treename="", TFile* file=0) 
		: SimulationOutput(),tree(0),dobunch(true) {
		
		// create a new ROOT file if necessary
		if(rootfile==0)
			if(file==0) rootfile = new TFile("merlin.root","RECREATE"); 
			else        rootfile = file;

		// create a tree 
		if(treename!="") NewTree(treename);

		// initialize blocks
		Output.AllFalse();
	}
	
        ~TrackingOutputROOT(){
		if(rootfile) rootfile->Write();
        }

	// Define which kind of variable block
	// block is one or many (additive) of the following
	// base         - z, reference momentum p0, #particles/slices in bunch
	// names        - type, name
	// bunchFirst   - means      - m_x,m_xp,m_y,m_yp,m_ct,m_dp   : <x>,<x'> etc.
	// bunchSecond  - sigmas     - s_x,s_xp,s_y,s_yp,s_ct,s_dp   : sigma(x) etc.
	// emittance    - emittances - gex,gey,geyc(energy corr.)
	// Twiss        - twiss parameters - ax,bx,dx,dxp,ay,by,dy,dyp
	// BPMCor       - BPM readings and Corrector settings - BPM_x BPM_y, XCor, YCor
	//                Only valid if BPMs etc. included e.g. in main to->AddIdentifier("BPM.*");
	// Cav          - Cavity phase and voltage
	//
	struct {
		bool base;
		bool names;
		bool bunchFirst;
		bool bunchSecond;
		bool emittance;
		bool Twiss;
		bool BPMCor;
		bool Cav;
		void AllFalse(){
			base=false;
			names=false;
			bunchFirst=false;
			bunchSecond=false;
			emittance=false;
			Twiss=false;
			BPMCor=false;
			Cav=false;
		};
	} Output;

	// Create a new TTree
	bool NewTree(const std::string& tname);

	// Access the current TTree directly
	TTree& GetTree() { return *tree; }

        // struct for transfering data
        typedef struct {

		Char_t type[16];// multiple of 4!
		Char_t name[16];// multiple of 4!

		Double_t z;
		Double_t p0;
		Int_t    n;

		Double_t m_x;
		Double_t m_y;
		Double_t m_xp;
		Double_t m_yp;
		Double_t m_dp;
		Double_t m_ct;

		Double_t gex;
		Double_t gey;
		Double_t geyc;

		Double_t ax;
		Double_t bx;
		Double_t ay;
		Double_t by;
		Double_t Dx;
		Double_t Dy;
		Double_t Dxp;
		Double_t Dyp;

		Double_t s_x;
		Double_t s_y;
		Double_t s_dp;
		Double_t s_ct;

		Double_t XCor;
		Double_t YCor;
		Double_t BPM_x;
		Double_t BPM_y;

		Double_t V;
		Double_t phi;

        } OutputStruct;

	void AddIdentifierAllMags();

	// dump bunch into root file at certain location
	// be aware that a location may not be unique (eg. BPM.*) -> many bunches 
	void DumpBunchAt(const string& ident){
		this->AddIdentifier(ident);
		dumpAt.push_back(ident);
	};
	int  BunchSize(const Bunch* bunch); // return the # of particles/slices of a bunch
	
	// switch bunch output on/off 
	void DoBunch(bool dob=true){dobunch=dob;};
	
	// additional utillities
	// 6d bunch struct for ROOT
        typedef struct bunch6D {
		Double_t x;
		Double_t xp;
		Double_t y;
		Double_t yp;
		Double_t ct;
		Double_t dp;
		Double_t Q;  // Q for SMP / Qtot/n for Particle 
        } OutputBunch;

	void DumpBunch(const Bunch* SB,const string& treename);

protected:

	virtual void Record(const ComponentFrame* frame, const Bunch* bunch);
	virtual void RecordInitialBunch(const Bunch* bunch) { zComponent=0; Record("INITIAL",bunch); }
        virtual void RecordFinalBunch(const Bunch* bunch){ Record("FINAL",bunch); }

private:
	void Record(const string&,const Bunch*);

	static TFile*  rootfile;     // the one and only ROOT file 
	TTree*             tree;     // and tree

	// for transfering data
	OutputStruct data;         
	OutputBunch  outputBunch;         


	double zComponent;

	// identifiers for bunch output
	vector<StringPattern> dumpAt;
	// true if in dumpAt
	bool isDumpAt(const string& id) {
		if(!dobunch) return false;
		for(vector<StringPattern>::const_iterator p=dumpAt.begin();p!=dumpAt.end();p++) {
			if(p->Match(id)) return true;
		}
		return false;
	};
	// fills dumpAt -  see usage
	void DumpBunchAt(const Bunch*, const std::string&);

	//map and counters for bunch names		
	map<string,int> bm;
	
	bool dobunch;

};

#endif

