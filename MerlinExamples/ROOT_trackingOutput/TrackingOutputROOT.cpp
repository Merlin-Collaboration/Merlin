#include "TrackingOutputROOT.h"
#include <fstream>
#include "BasicTransport/NormalTransform.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include <iomanip>

#include "BasicTransport/RMap.h"
#include "AcceleratorModel/ActiveMonitors/BPM.h"
#include "AcceleratorModel/Components.h"

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;

namespace {

class BPMBuffer : public BPM::Buffer {
public:
   BPMBuffer(double*ct,double*x,double*y) : CT(ct),X(x),Y(y) {}
   void Record(const BPM& theBPM, const BPM::Data& theData);
private:
   double* X;
   double* Y;
   double* CT;
};
void BPMBuffer::Record(const BPM& theBPM, const BPM::Data& theData){
   *X=theData.x.value;
   *Y=theData.y.value;
   *CT=theData.ct;
}

void SigmaMatrixToBeamData0(const PSmoments& S0, BeamData& t);

double DispersionCorrectedEmittance(const PSmoments& S,size_t a,size_t ap)	{
		double s36 = S(a,ps_DP);
		double s46 = S(ap,ps_DP);
		double dp2 = S.var(ps_DP);
		double s33 = S.var(a)-s36*s36/dp2;
		double s34 = S(a,ap)-s36*s46/dp2;
		double s44 = S.var(ap)-s46*s46/dp2;
		return sqrt(s33*s44-s34*s34);
};

}
void TrackingOutputROOT::Record(const ComponentFrame* frame, const Bunch* bunch)
{
    if(frame->IsComponent()) {
	zComponent=frame->GetPosition();
	string id = (*frame).GetComponent().GetQualifiedName();

	data.YCor=data.BPM_y=0;
	data.XCor=data.BPM_x=0;
	data.phi=data.V=0;

	//BPMs
	if(id.find("BPM.")==0) {
		const BPM* bb = static_cast<const  BPM* >(&(frame->GetComponent()));
		BPM b =*bb;
		double ct;
		BPM::Buffer* newBuf = new BPMBuffer(&ct,&data.BPM_x,&data.BPM_y); // delete?
		b.AddBuffer(newBuf);
		b.MakeMeasurement(*bunch);
		//cout<<ct<<" "<<data.bpm_Y<<endl;
	}
	if(id.find("YCor.")==0) {
		const YCor* ycor = static_cast<const  YCor* >(&(frame->GetComponent()));
		//cout<<"YCor "<<ycor->GetFieldStrength()<<endl;
		data.YCor=ycor->GetFieldStrength();
	}
	if(id.find("XCor.")==0) {
		const XCor* xcor = static_cast<const  XCor* >(&(frame->GetComponent()));
		//cout<<"XCor "<<xcor->GetFieldStrength()<<endl;
		data.XCor=xcor->GetFieldStrength();
	}
	if(id.find("TWRFStructure.")==0) {
		const TWRFStructure* lcav = static_cast<const  TWRFStructure* >(&(frame->GetComponent()));
		data.phi=lcav->GetPhase();
		data.V=lcav->GetVoltage();
	}
	Record(id,bunch);
    }
}

void TrackingOutputROOT::Record(const string& id, const Bunch* bunch) {

	PSmoments S;
	bunch->GetMoments(S);
	double E = bunch->GetReferenceMomentum();
	       E*=1+S.mean(ps_DP);

	double sigX  = sqrt(S.var(ps_X));
	double sigY  = sqrt(S.var(ps_Y));
	double sigDP = sqrt(S.var(ps_DP));
	double sigCT = sqrt(S.var(ps_CT));
        
        BeamData  B;
        SigmaMatrixToBeamData0(S,B);

	double gamma = E/MeV/ElectronMassMeV;
	double gex   = gamma*ProjectedEmittance(S,ps_X,ps_XP);
	double gey   = gamma*ProjectedEmittance(S,ps_Y,ps_YP);
	double gexc  = S.var(ps_DP) != 0 ? gamma*DispersionCorrectedEmittance(S,ps_X,ps_XP) : gex;
	double geyc  = S.var(ps_DP) != 0 ? gamma*DispersionCorrectedEmittance(S,ps_Y,ps_YP) : gey;
 
	size_t n = id.find('.');
        strncpy(data.type,id.substr(0,n).c_str(),15);data.type[15]='\0';
        strncpy(data.name,id.substr(n+1).c_str(),15);data.type[15]='\0';

//	data.z     = bunch->GetReferenceTime(); //S - nout unique starts at 0 for each tracker
	data.z     = zComponent;
        data.p0    = bunch->GetReferenceMomentum();;
	data.n     = BunchSize(bunch);
	data.m_x   = B.x0;
	data.m_y   = B.y0;

	data.m_dp   = B.p0; // that's dp_0 not the reference momentum
	data.m_ct   = B.ct0;

	data.m_xp  = B.xp0;
	data.m_yp  = B.yp0;
	data.gex   = gex;
	data.gey   = gey;
	data.geyc  = geyc;
	
	data.ax  = B.alpha_x;
	data.bx  = B.beta_x;
	data.ay  = B.alpha_y;
	data.by  = B.beta_y;

	data.Dx  = B.Dx;
	data.Dy  = B.Dy;
	data.Dxp  = B.Dxp;
	data.Dyp  = B.Dyp;

	data.s_x   = sigX; // sqrt(var)
	data.s_y   = sigY;
	data.s_dp  = sigDP;
	data.s_ct  = sigCT;
	
	if(tree) tree->Fill(); // data -> ROOT tree
	else {
		cout<<"TrackingOutputROOT: No root tree defined! "<<endl;
		abort();
	}
	if( isDumpAt(id) ) DumpBunchAt(bunch,id);
}

bool TrackingOutputROOT::NewTree(const std::string& tname) {
	
	if(tree!=0){ 
                tree->Write();
                delete tree;
                bm.clear();
        }

	tree = new TTree(tname.c_str(),"TrackingTree");
        if(tree) {
	   int Bufsize=100000;
	   
	   // block names
	   if(Output.names){
	        tree->Branch("type", &data.type, "type[16]/C",Bufsize);
        	tree->Branch("name", &data.name, "name[16]/C",Bufsize);
           }
           // block base
	   if(Output.base){
		tree->Branch("z",    &data.z,    "z/D",Bufsize);
		tree->Branch("p0",   &data.p0,   "p0/D",Bufsize);
		tree->Branch("n",    &data.n,    "n/I",Bufsize);
	   }
	   // block bunchFirst
	   if(Output.bunchFirst){
		tree->Branch("m_x",  &data.m_x,  "m_x/D",Bufsize);
		tree->Branch("m_y",  &data.m_y,  "m_y/D",Bufsize);
		tree->Branch("m_xp", &data.m_xp, "m_xp/D",Bufsize);
		tree->Branch("m_yp", &data.m_yp, "m_yp/D",Bufsize);
		tree->Branch("m_dp", &data.m_dp, "m_dp/D",Bufsize);
		tree->Branch("m_ct", &data.m_ct, "m_ct/D",Bufsize);
	   }
	   // block emittance
	   if(Output.emittance){
		tree->Branch("gex",  &data.gex,  "gex/D",Bufsize);
		tree->Branch("gey",  &data.gey,  "gey/D",Bufsize);
		tree->Branch("geyc", &data.geyc, "geyc/D",Bufsize);
	   }
	   // block Twiss
	   if(Output.Twiss){
		tree->Branch("ax",   &data.ax,   "ax/D",Bufsize);
		tree->Branch("bx",   &data.bx,   "bx/D",Bufsize);
		tree->Branch("ay",   &data.ay,   "ay/D",Bufsize);
		tree->Branch("by",   &data.by,   "by/D",Bufsize);
		tree->Branch("Dx",   &data.Dx,   "Dx/D",Bufsize);
		tree->Branch("Dy",   &data.Dy,   "Dy/D",Bufsize);
		tree->Branch("Dxp",   &data.Dxp,   "Dxp/D",Bufsize);
		tree->Branch("Dyp",   &data.Dyp,   "Dyp/D",Bufsize);
	   }
	   // block bunchSecond
	   if(Output.bunchSecond){
		tree->Branch("s_x", &data.s_x, "s_x/D",Bufsize);
		tree->Branch("s_y", &data.s_y, "s_y/D",Bufsize);
		tree->Branch("s_dp",&data.s_dp,"s_dp/D",Bufsize);
		tree->Branch("s_ct",&data.s_ct,"s_ct/D",Bufsize);
	   }
	   // block BPMCor
	   if(Output.BPMCor){
		tree->Branch("BPM_x",&data.BPM_x,"BPM_x/D",Bufsize);
		tree->Branch("BPM_y",&data.BPM_y,"BPM_y/D",Bufsize);
		tree->Branch("XCor", &data.XCor, "XCor/D",Bufsize);
		tree->Branch("YCor", &data.YCor, "YCor/D",Bufsize);
	   }
	   // block Cav
	   if(Output.names){
		tree->Branch("V", &data.V, "V/D",Bufsize);
		tree->Branch("phi", &data.phi, "phi/D",Bufsize);
	   }
	   
           return true;
	}
	                                               
        return false;
};

namespace {
void SigmaMatrixToBeamData0(const PSmoments& S0, BeamData& t) {

    PSmoments S=S0;
    t=BeamData();

    // First check if we have dispersion:
    double d2 = S.var(ps_DP);
    if(d2!=0) {
        t.Dx = S(ps_X,ps_DP)/d2;
        t.Dxp = S(ps_XP,ps_DP)/d2;
        t.Dy = S(ps_Y,ps_DP)/d2;
        t.Dyp = S(ps_YP,ps_DP)/d2;

        // Remove energy correlations
        RMap D(DispersionMatrix(-t.Dx,-t.Dxp,-t.Dy,-t.Dyp));
        D.Apply(S);
    }

    t.emit_x = ProjectedEmittance(S,ps_X,ps_XP);
    t.emit_y = ProjectedEmittance(S,ps_Y,ps_YP);

    // Now calculate beta, alpha,
    t.beta_x = S.var(ps_X)/t.emit_x;
    t.beta_y = S.var(ps_Y)/t.emit_y;
    t.alpha_x = -S(ps_X,ps_XP)/t.emit_x;
    t.alpha_y = -S(ps_Y,ps_YP)/t.emit_y;

    // Remove the correlation
    //RealMatrix B=InverseBetaTransform(t.beta_x,t.beta_y,t.alpha_x,t.alpha_y);
    //RMap(B).Apply(S);

    t.sig_z  = S0.std(ps_CT);
    t.sig_dp = S0.std(ps_DP);

    // Centroid
    t.x0  = S0[0];
    t.xp0 = S0[1];
    t.y0  = S0[2];
    t.yp0 = S0[3];
    t.ct0 = S0[4];
    t.p0  = S0[5]; // not actually the energy, but dp/p!

}

}

TFile* TrackingOutputROOT::rootfile = 0;

void TrackingOutputROOT::AddIdentifierAllMags(){

	AddIdentifier("SectorBend.*");
	AddIdentifier("Quadrupole.*");
	AddIdentifier("SkewQuadrupole.*");
	AddIdentifier("Sextupole.*");
	AddIdentifier("SkewSextupole.*");
	AddIdentifier("Octupole.*");
	AddIdentifier("Decapole.*");
	AddIdentifier("Solenoid.*");

};

// bunch utillities
void TrackingOutputROOT::DumpBunchAt(const Bunch* bunch,const string& id){

	// add current tree name + name part of id 
	// (+ counter if id appears more than onece)
	stringstream btreename;
	btreename<<tree->GetName()<<"_";
	string s=id.substr(id.rfind(".")+1);
	if(bm[s]++==0) btreename<<s; 
	else           btreename<<s<<"_"<<bm[s]-1;

	DumpBunch(bunch,btreename.str());
}
int TrackingOutputROOT::BunchSize(const Bunch* bunch){
	if (  typeid(*bunch)==typeid(SMPBunch) ) {
		const SMPBunch* SB = static_cast<const SMPBunch*>(bunch);
		return SB->Size();
	} else {
		const ParticleBunch* PB = static_cast<const ParticleBunch*>(bunch);
		return PB->size();
	}
};
// write bunch tree as - can also be called externally
void TrackingOutputROOT::DumpBunch(const Bunch* bunch,const string& btreename){
	if(bunch==0) return;
	if(btreename=="") {
		cout<<"TrackingOutputROOTROOT::DumpBunch: without name ?"<<endl;
		return;
	}
	if (  typeid(*bunch)==typeid(SMPBunch) ) {
		const SMPBunch* SB = static_cast<const SMPBunch*>(bunch);
		TTree tree(btreename.c_str(),"SMPBunchTree");
		tree.Branch("Bunch",&outputBunch,"x/D:xp/D:y/D:yp/D:ct/D:dp/D:Q/D");
		for(SMPBunch::const_iterator sp=SB->begin(); sp!=SB->end(); sp++) {
			outputBunch.x  = sp->x();
			outputBunch.xp = sp->xp();
			outputBunch.y  = sp->y();
			outputBunch.yp = sp->yp();
			outputBunch.ct = sp->ct();
			outputBunch.dp = sp->dp();
			outputBunch.Q  = sp->Q();
			tree.Fill();
            	}
            	tree.Write();
	} else {
		const ParticleBunch* PB = static_cast<const ParticleBunch*>(bunch);
		double q=PB->GetTotalCharge()/PB->size();
		TTree tree(btreename.c_str(),"ParticleBunchTree");
		tree.Branch("Bunch",&outputBunch,"x/D:xp/D:y/D:yp/D:ct/D:dp/D:Q/D");
		for(ParticleBunch::const_iterator pb = PB->begin(); pb!= PB->end(); pb++){
			outputBunch.x  = pb->x();
			outputBunch.xp = pb->xp();
			outputBunch.y  = pb->y();
			outputBunch.yp = pb->yp();
			outputBunch.ct = pb->ct();
			outputBunch.dp = pb->dp();
			outputBunch.Q  = q;
			tree.Fill();
            	}
            	tree.Write();
       }
};
