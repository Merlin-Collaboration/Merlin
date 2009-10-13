/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/09/15 13:43:32 $
// $Revision: 1.11 $
// 
/////////////////////////////////////////////////////////////////////////

#include "MADInterface/XTFFInterface.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cctype>

#include "AcceleratorModel/Components.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/Supports/SupportStructure.h"
#include "AcceleratorModel/Supports/MagnetMover.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "MADInterface/ConstructSrot.h"
#include "Exception/MerlinException.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;
using namespace std;

// data structure for XTFF data
struct XTFFInterface::XTFF_Data {
	std::string keywrd;
	std::string label;
	std::string type;
	std::string fdn;
	double realdat[11];
	double operator[](size_t n) const { return realdat[n];}
};

namespace {

	typedef XTFFInterface::XTFF_Data Data;

	double energy;      // current energy
	double beamload;        // energy loss due to beamloading

	double z_total; // current total length
	double Qt;      // charge for beamloading

	// skip n input lines
	void SkipLines(istream& is, int n)
	{
		string dummy;
		while(is && n--)
			getline(is,dummy);
	}

	// return real value in columnds c1-c2
	double RealValue(const string& dat, int c1, int c2)
	{
		return atof(dat.substr(c1-1,c2-c1+1).c_str());
	}

	// return string value in cols c1-c2
	string StringValue(const string& dat, int c1, int c2)
	{
		string rv = dat.substr(c1-1,c2-c1+1);
		int n = rv.find_first_of(' ');
		return n==string::npos ? rv : rv.substr(0,n);
	}


	// Parse single element record
	void ParseXTFF(istream& is, Data& data)
	{
		string ipline;

		// first line
		getline(is,ipline);
		data.keywrd = StringValue(ipline,1,4);
		data.label = StringValue(ipline,5,20);
		data.type = StringValue(ipline,98,113);

		data.realdat[0] = RealValue(ipline,21,32);
		data.realdat[1] = RealValue(ipline,33,48);
		data.realdat[2] = RealValue(ipline,49,64);
		data.realdat[3] = RealValue(ipline,65,80);
		data.realdat[4] = RealValue(ipline,81,96);
		data.realdat[5] = RealValue(ipline,115,130);

		// second line
		getline(is,ipline);

		data.realdat[6] = RealValue(ipline,1,16);
		data.realdat[7] = RealValue(ipline,17,32);
		data.realdat[8] = RealValue(ipline,33,48);
		data.realdat[9] = RealValue(ipline,49,64);
		data.realdat[10] = RealValue(ipline,65,80);

		if(ipline.length()>82)
			data.fdn = StringValue(ipline,82,105);
	}

	// parameter keyword array locations

#define L      0
#define ANGLE  1
#define K1     2
#define K2     3
#define K3     7
#define K4     7
#define ENERGY 5
#define APER   4
#define HGAP   4
#define TILT   6
#define FREQ   7
#define VOLT   8
#define LAG    9
#define ELOSS 10
#define KS     7
#define E1     7
#define E2     8
#define H1     9
#define H2    10
#define HKICK  6
#define VKICK  7
#define XGAP   6
#define YGAP   7
#define SROT   7

	Drift* ConstructDrift(const Data& data)
	{
		return new Drift(data.label,data[L]);
	}

	SectorBend* ConstructSectorBend(const Data& data)
	{
		double len   = data[L];
		double angle = data[ANGLE];
		double k1    = data[K1];
		double k2    = data[K2];
		double e1    = data[E1];
		double e2    = data[E2];
		double tilt  = data[TILT];
		double brho  = energy/eV/SpeedOfLight;
		double hg    = data[HGAP];
		double h = angle/len;

		SectorBend* bend = new SectorBend(data.label,len,h,brho*h);

		if(k1!=0)  // mixed function dipole
			bend->SetB1(brho*k1);

		if(tilt!=0)
			(*bend).GetGeometry().SetTilt(tilt);

		// add pole-face rotation information
		SectorBend::PoleFace* entrPF = (e1!=0)? new SectorBend::PoleFace(e1,0,hg) : 0;
		SectorBend::PoleFace* exitPF = (e2!=0)? new SectorBend::PoleFace(e2,0,hg) : 0;
		bend->SetPoleFaceInfo(entrPF,exitPF);

		return bend;
	}

	SectorBend* ConstructRectBend(const Data& data)
	{
		double len   = data[L];
		double angle = data[ANGLE];
		double k1    = data[K1];
		double k2    = data[K2];
		double e1    = angle/2;
		double e2    = angle/2;
		double tilt  = data[TILT];
		double brho  = energy/eV/SpeedOfLight;
		double hg    = data[HGAP];

		double h=0;
		if(angle!=0) {
			h = 2*sin(angle/2)/len;
			len = angle/h;
		}


		SectorBend* bend = new SectorBend(data.label,len,h,brho*h);

		if(k1!=0)  // mixed function dipole
			bend->SetB1(brho*k1);

		if(tilt!=0)
			(*bend).GetGeometry().SetTilt(tilt);

		// add pole-face rotation information
		SectorBend::PoleFace* entrPF = (e1!=0)? new SectorBend::PoleFace(e1,0,hg) : 0;
		SectorBend::PoleFace* exitPF = (e2!=0)? new SectorBend::PoleFace(e2,0,hg) : 0;
		bend->SetPoleFaceInfo(entrPF,exitPF);

		return bend;
	}

	Quadrupole* ConstructQuadrupole(const Data& data)
	{
		double len   = data[L];
		double k1    = data[K1];
		double tilt  = data[TILT];
		assert(tilt==0);
		double brho  = energy/eV/SpeedOfLight;
		return new Quadrupole(data.label,len,brho*k1);
	}

	SkewQuadrupole* ConstructSkewQuadrupole(const Data& data)
	{
		double len   = data[L];
		double k1    = data[K1];
		double tilt  = data[TILT];
		double brho  = energy/eV/SpeedOfLight;
		return new SkewQuadrupole(data.label,len,brho*k1);
	}

	Sextupole* ConstructSextupole(const Data& data)
	{
		double len   = data[L];
		double k2    = data[K2];
		double tilt  = data[TILT];
		assert(tilt==0);
		double brho  = energy/eV/SpeedOfLight;
		return new Sextupole(data.label,len,brho*k2);
	}

	SkewSextupole* ConstructSkewSextupole(const Data& data)
	{
		double len   = data[L];
		double k2    = data[K2];
		double tilt  = data[TILT];
		double brho  = energy/eV/SpeedOfLight;
		return new SkewSextupole(data.label,len,brho*k2);
	}

	Octupole* ConstructOctupole(const Data& data)
	{
		double len   = data[L];
		double k3    = data[K3];
		double tilt  = data[TILT];
		assert(tilt==0);
		double brho  = energy/eV/SpeedOfLight;
		return new Octupole(data.label,len,brho*k3);
	}

	Decapole* ConstructDecapole(const Data& data)
	{
		double len   = data[L];
		double k4    = data[K4];
		double tilt  = data[TILT];
		assert(tilt==0);
		double brho  = energy/eV/SpeedOfLight;
		return new Decapole(data.label,len,brho*k4);
	}


	XCor* ConstructXCor(const Data& data)
	{
		return new XCor(data.label,data[L]);
	}

	YCor* ConstructYCor(const Data& data)
	{
		return new YCor(data.label,data[L]);
	}

	TWRFStructure* ConstructCavity(const Data& data)
	{
		double len = data[L];
		double volt = data[VOLT]*MV;
		double phase = twoPi*data[LAG];
		double freq = data[FREQ]*MHz;
		double eloss = data[ELOSS]*Volt;

		// update energy
		double bloading = eloss*Qt;
		double dE = volt*cos(phase);
		/*
		using std::setw;
		cout<<setw(12)<<data.label.c_str();
		cout<<setw(10)<<fixed<<setprecision(3)<<dE/MV;
		cout<<setw(10)<<fixed<<setprecision(3)<<bloading/MV<<endl;
		*/
		energy += (dE-bloading)/GeV;
		beamload += bloading/GeV;
		return new TWRFStructure(data.label,len,freq,volt/len,phase);
	}

	TransverseRFStructure* ConstructCrabCavity(const Data& data, char c)
	{
		double len = data[L];
		double volt = data[VOLT]*MV;
		double phase = twoPi*data[LAG];
		double freq = data[FREQ]*MHz;
		double roll;

		switch(c) {
case 'H': roll = 0; break;
case 'V': roll = pi/2.0; break;
default: 
	throw MerlinException("Cannot deduce orientation for crab cavity");
	break;
		}
		return new TransverseRFStructure(data.label,len,freq,volt/len,phase,roll);
	}

	BPM* ConstructBPM(const Data& data)
	{
		return new BPM(data.label,data[L]);
	}

	RMSProfileMonitor* ConstructProfileMonitor(const Data& data)
	{
		return new RMSProfileMonitor(data.label,data[L]);
	}

	Solenoid* ConstructSolenoid(const Data& data)
	{
		double brho = energy/eV/SpeedOfLight;
		return new Solenoid(data.label,data[L],brho*data[KS]);
	}

	Marker* ConstructMarker(const Data& data)
	{
		return new Marker(data.label);
	}

#define SKQ_TILT 0.78539816
#define SKS_TILT 0.52359878

	inline bool is_skewquad(double tilt) { return fabs(tilt/SKQ_TILT-1.0)<1e-03;}
	inline bool is_skewsext(double tilt) { return fabs(tilt/SKS_TILT-1.0)<1e-03;}


#define TYPEIS(kw) (dat.keywrd == #kw)

}; // end namespace

void XTFFInterface::ConstructComponent(XTFF_Data& dat)
{
	if(TYPEIS(HMON)||TYPEIS(VMON))
		dat.keywrd = "MONI";
	//      else if(TYPEIS(MARK))
	//              dat.keywrd = "DRIF";

	if(driftTypes.find(dat.keywrd)!=driftTypes.end()) {
		cerr<<"WARNING: treating "<<dat.keywrd<<" as DRIFT"<<endl;
		dat.keywrd = "DRIF";
	}

	AcceleratorComponent* c;

	if(dat.keywrd=="DRIF")
		c = mc->AppendComponent(ConstructDrift(dat));
	else if(dat.keywrd=="QUAD") {
		if(is_skewquad(dat[TILT]))
			c = mc->AppendComponent(ConstructSkewQuadrupole(dat));
		else
			c = mc->AppendComponent(ConstructQuadrupole(dat));
	}
	else if(dat.keywrd=="SBEN")
		c = mc->AppendComponent(ConstructSectorBend(dat));
	else if(dat.keywrd=="RBEN")
		c = mc->AppendComponent(ConstructRectBend(dat));
	else if(dat.keywrd=="SEXT") {
		if(is_skewsext(dat[TILT]))
			c = mc->AppendComponent(ConstructSkewSextupole(dat));
		else
			c = mc->AppendComponent(ConstructSextupole(dat));
	}
	else if(dat.keywrd=="OCTU")
		c = mc->AppendComponent(ConstructOctupole(dat));
	else if(dat.keywrd=="DECA")
		c = mc->AppendComponent(ConstructDecapole(dat));
	else if(dat.keywrd=="LCAV") {
		if(dat.label.substr(0,4)=="CRAB")
			c = mc->AppendComponent(ConstructCrabCavity(dat,(dat.label)[4]));
		else
			c = mc->AppendComponent(ConstructCavity(dat));
	}
	else if(dat.keywrd=="SOLE")
		c = mc->AppendComponent(ConstructSolenoid(dat));
	else if(dat.keywrd=="HKIC")
		c = mc->AppendComponent(ConstructXCor(dat));
	else if(dat.keywrd=="VKIC")
		c = mc->AppendComponent(ConstructYCor(dat));
	else if(dat.keywrd=="MONI")
		c = mc->AppendComponent(ConstructBPM(dat));
	else if(dat.keywrd=="WIRE")
		c = mc->AppendComponent(ConstructProfileMonitor(dat));
	else if(dat.keywrd=="MARK") {
		if(girders && dat.label[1]=='_') {
			string girderName = dat.label;
			if(!frameStack.empty() && girderName == frameStack.top()) {
				if(logos)
					(*logos)<<"<--- "<<girderName<<" END "<<mc->GetCurrentFrameDepth()<<endl;
				mc->EndFrame();
				frameStack.pop();
			}
			else {
				frameStack.push(girderName);
				ConstructNewFrame(girderName);
			}
			c=0;
		}
		else if(dat.label=="VPIV") {
			// here we inserted a vertical kink in the beamline
			// which we will later use to allow terrain following
			// We also include a marker for reference.
			mc->AppendComponentFrame(ConstructXrot(0.0,"VPIV"));
			c = mc->AppendComponent(ConstructMarker(dat));
		}
		else
			c = mc->AppendComponent(ConstructMarker(dat));
	}
	else if(dat.keywrd=="RCOL") {
		// construct a drift with a rectangular aperture
		c = mc->AppendComponent(ConstructDrift(dat));
		if(incApertures)
			c->SetAperture(new RectangularAperture(2*dat[XGAP],2*dat[YGAP]));
	}
	else if(dat.keywrd=="SROT") {
		mc->AppendComponentFrame(ConstructSrot(dat[SROT],dat.label));
		c=0;
	}
	else {
		cerr<<"WARNING: treating "<<dat.keywrd<<" as DRIFT"<<endl;
		c = mc->AppendComponent(ConstructDrift(dat));
	}

	if(c && incApertures && dat[APER]!=0) {
		c->SetAperture(new CircularAperture(dat[APER]));
	}

	if(c) {
		z_total += c->GetLength();
		if(logos) {
			(*logos)<<setw(10)<<left<<(*c).GetName().c_str();
			(*logos)<<setw(16)<<left<<(*c).GetType().c_str();
			(*logos)<<setw(10)<<right<<fixed<<setprecision(3)<<z_total;
			(*logos)<<setw(10)<<right<<fixed<<setprecision(3)<<energy;
			(*logos)<<setw(10)<<right<<fixed<<setprecision(3)<<beamload;
			(*logos)<<setw(10)<<right<<fixed<<setprecision(3)<<dat[ENERGY];
			(*logos)<<endl;
		}
	}
}

void XTFFInterface::ConstructNewFrame (const string& name)
{
	SequenceFrame* newFrame;

	switch( name[0] ) {
		case 'F':
			newFrame = new SequenceFrame(name.substr(2));
			break;
		case 'S':
			newFrame = new SimpleMount(name.substr(2));
			break;
		case 'G':
			// newFrame = new SimpleMount(name.substr(2));
			newFrame = new GirderMount(name.substr(2));
			break;
		case 'M':
			newFrame = new MagnetMover(name.substr(2));
			break;
		default:
			throw MerlinException(string("Unknown frame key: ")+name);
			break;
	}

	mc->NewFrame(newFrame);
	if(logos)
		(*logos)<<"---> "<<name<<" BEGIN "<<mc->GetCurrentFrameDepth()<<endl;
}

XTFFInterface::XTFFInterface(const string& fname, double Nb, ostream* logstream)
: ifs(new ifstream(fname.c_str())),mc(0),beam0(0),nb(Nb),logos(logstream),incApertures(true),
girders(false),frameStack()
{
	if(!(*ifs)) {
		string msg = "cannot open file "+fname;
		delete ifs;
		throw runtime_error(msg);
	}
}

XTFFInterface::XTFFInterface(double Nb, ostream* logstream)
: ifs(0),mc(0),beam0(0),nb(Nb),logos(logstream),incApertures(true),
girders(false),frameStack()
{}


XTFFInterface::~XTFFInterface()
{
	if(ifs) delete ifs;
	if(mc!=0) delete mc;
	if(beam0!=0) delete beam0;
}

pair<AcceleratorModel*,BeamData*> XTFFInterface::Parse()
{
	if(beam0!=0)
		delete beam0;
	beam0 = new BeamData();
	int nelm = ParseHeader(beam0);

	Parse1(nelm,true);
	return GetModel();
}

pair<AcceleratorModel*,BeamData*> XTFFInterface::Parse(double P_ref)
{
	if(beam0!=0)
		delete beam0;
	beam0 = new BeamData();
	int nelm = ParseHeader(beam0);
	energy = P_ref;
	Parse1(nelm,true);
	return GetModel();
}

void XTFFInterface::AppendModel(const string& fname)
{
	if(ifs)
		delete ifs;
	ifs = new ifstream(fname.c_str());

	if(!(*ifs)) {
		string msg = "cannot open file "+fname;
		delete ifs;
		throw runtime_error(msg);
	}

	BeamData* b = new BeamData();
	int nelm = ParseHeader(b);

	if(beam0==0)
		beam0 = b;
	else
		delete b;

	Parse1(nelm,false);
}


void XTFFInterface::Parse1(int nelm, bool isNewModel)
{
	if(mc==0 || isNewModel) {
		if(mc!=0)
			delete mc;
		mc = new AcceleratorModelConstructor();
		mc->NewModel();
		z_total=0;
		beamload=0;
		energy = beam0->p0;
	}

	if(logos)
		(*logos)<<"Initial beam energy: "<<energy<<" GeV"<<endl;

	while((*ifs) && nelm--) {
		XTFF_Data dat;
		ParseXTFF((*ifs),dat);
		ConstructComponent(dat);
		SkipLines((*ifs),3);
	}
}

pair<AcceleratorModel*,BeamData*> XTFFInterface::GetModel()
{
	if(logos) {
		(*logos)<<endl;
		mc->ReportStatistics(*logos);
		(*logos)<<"\nFinal ";
		if(Qt!=0)
			(*logos)<<"(loaded) ";
		(*logos)<<"beam energy: "<<energy<<" GeV"<<endl;
	}

	BeamData* beam = beam0;
	beam0=0;
	return make_pair(mc->GetModel(),beam);
}

int XTFFInterface::ParseHeader(BeamData* bdat)
{
	string ipline;

	// extract number of elements from first line
	getline((*ifs),ipline);
	int n = RealValue(ipline,57,65);

	// Skip next header line
	getline((*ifs),ipline);

	XTFF_Data dat;
	ParseXTFF((*ifs),dat); // initial element
	bdat->p0 = dat[ENERGY];
	bdat->charge = nb;
	Qt = (bdat->charge)*ElectronCharge;

	// Now extract initial beam conditions from
	// next three records

	double dummy; // used for mux,muy

	(*ifs)>>bdat->alpha_x;
	(*ifs)>>bdat->beta_x;
	(*ifs)>>dummy; // mux
	(*ifs)>>bdat->Dx;
	(*ifs)>>bdat->Dxp;

	(*ifs)>>bdat->alpha_y;
	(*ifs)>>bdat->beta_y;
	(*ifs)>>dummy; // muy
	(*ifs)>>bdat->Dy;
	(*ifs)>>bdat->Dyp;

	(*ifs)>>bdat->x0;
	(*ifs)>>bdat->xp0;
	(*ifs)>>bdat->y0;
	(*ifs)>>bdat->yp0;

	// remove remainder of this line (suml)
	getline((*ifs),ipline);

	return n-1;
}

void XTFFInterface::TreatTypeAsDrift(const string& dt)
{
	char buff[5];
	for(size_t i=0; i<4; i++)
		buff[i]=toupper(dt[i]);
	buff[4]=0;
	driftTypes.insert(buff);
}
