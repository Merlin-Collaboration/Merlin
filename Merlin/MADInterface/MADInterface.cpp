/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/01/31 16:20:21 $
// $Revision: 1.15 $
//
/////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include "AcceleratorModel/Components.h"
#include "AcceleratorModel/StdComponent/Spoiler.h"

#include "IO/MerlinIO.h"
#include "NumericalUtils/utils.h"

// SupportStructure
#include "AcceleratorModel/Supports/SupportStructure.h"
// MagnetMover
#include "AcceleratorModel/Supports/MagnetMover.h"
// SequenceFrame
#include "AcceleratorModel/Frames/SequenceFrame.h"
// SimpleApertures
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/Apertures/TiltedAperture.hpp"

// AcceleratorModelConstructor
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
// PhysicalConstants
#include "NumericalUtils/PhysicalConstants.h"
// MADKeyMap
#include "MADInterface/MADKeyMap.h"
// MADInterface
#include "MADInterface/MADInterface.h"
#include "MADInterface/ConstructSrot.h"
#include "Collimators/ResistiveWakePotentials.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

const char* material_names[] = {"Be", "C", "Al", "Cu", "W", "Pb"};

namespace {
// stack used to match MAD LINE pairs
// Note: we need this because we no longer preserve the complete MAD beamline structure.
stack<string> frameStack;

// function to calculate energy loss in dipole due to SR (returns dE in GeV)
inline double SRdE(double h, double len, double E)
{
	static const double Cg = 8.85e-05/twoPi;
	return Cg*pow(E,4)*h*h*len;
}

// utility routines for reading in MAD optics file
void Log(const string& tag, int depth, ostream& os)
{
	static const char* tab = "----|";
	while(depth--) os<<tab;
	os<<' '<<tag<<endl;
}

void StripHeader(istream& is)
{
	char buffer[512];
	char c;
	while(true && is)
	{
		is.get(c);
		if(c=='*'||c=='$'||c=='@')
		{
			is.getline(buffer,512);
		}
		else
		{
			is.putback(c);
			break;
		}
	}
}

inline string StripQuotes(const string& text)
{
	return text.substr(1,text.length()-2);
}




Aperture* ConstructAperture(const string& apstr)
{
	Aperture* ap;
	double w,h,d,t;
	string::size_type n,n2;

	switch(apstr[0])
	{
	case 'D':
		d = atof(apstr.substr(1).c_str())*millimeter;
		ap = new CircularAperture(d/2);
		break;

	case 'X':
		n = apstr.find_first_of('Y',1);
		assert(n!=string::npos);
		w = atof(apstr.substr(1,n-1).c_str())*millimeter;
		h = atof(apstr.substr(n+1).c_str())*millimeter;
		ap = new RectangularAperture(w,h);
		break;

	case 'T':
		n = apstr.find_first_of('X',1);
		n2 = apstr.find_first_of('Y',1);
		assert(n!=string::npos);
		t = atof(apstr.substr(1,n-1).c_str());
		w = atof(apstr.substr(n+1,n2-1).c_str())*millimeter;
		h = atof(apstr.substr(n2+1).c_str())*millimeter;
		// cout << "found a tilted aperture w="<<w<<" h="<<h<<" t="<<t<<endl;
		ap = new TiltedAperture(w,h,t);
		break;

		case '~':
		        ap = 0;
		        break;

		default:
			MERLIN_ERR<<"WARNING: unknown aperture definition ("<<apstr<<") ignored"<<endl;
			ap=0;
	};

	return ap;
}



void check_column_heading(istream& is, const string& hd)
{
	string s;
	is>>s;
	if(s!=hd)
	{
		MerlinIO::error() << "Bad file format: expected " << hd << " found " << s << endl;
		abort();
	}
}

}; // Namespace end

// Class MADInterface
MADInterface::MADInterface (const std::string& madFileName, double P0)
        : energy(P0),ifs(madFileName.empty() ? 0 : new ifstream(madFileName.c_str())),
        log(MerlinIO::std_out),logFlag(false),flatLattice(false),honMadStructs(false),
        incApertures(false),inc_sr(false),ctor(0),prmMap(0)
{
	if(ifs)
	{
		if(!(*ifs))
		{
			throw MerlinException(string("ERROR opening file ")+string(madFileName));
		}

		Initialise();
	}

	// By default, we currently treat the following MAD types as drifts
	TreatTypeAsDrift("MARKER"); // merlin bug!
	TreatTypeAsDrift("INSTRUMENT"); // merlin bug!
}

void MADInterface::Initialise()
{
	//The following checks were for a file that is not valid TFS format, changing to check for what MADX actually produces - JM
	//check_column_heading((*ifs),"*");
	//check_column_heading((*ifs),"NAME");
	//check_column_heading((*ifs),"KEYWORD");
	//Is MAD8 different?

	//check_column_heading((*ifs),"@");
	//check_column_heading((*ifs),"NAME");

	if(prmMap!=0)
	{
		delete prmMap;
	}

	string s;
	char c;
	bool tfs = 0;
	//This just grabs the top line, we really want the line that starts with *, aka the one with the headers.
	//getline((*ifs),s);
	//Fixed - JM:

	while((*ifs).good())
	{
		c = (*ifs).get();
		if(c=='*')
		{
			check_column_heading((*ifs),"NAME");
			check_column_heading((*ifs),"KEYWORD");
			tfs = 1;
			getline((*ifs),s);
			break;
		}
	}
	if(tfs == 0)
	{
		cout << "No Suitable TFS file headers found" << endl;
		abort();
	}

	prmMap = new MADKeyMap(s);
}

void MADInterface::AppendModel (const string& fname, double Pref)
{
    if(ifs)
        delete ifs;

    ifs =  new ifstream(fname.c_str());
    if(!(*ifs)) {
        MERLIN_ERR<<"ERROR opening file "<<fname<<endl;
        delete ifs;
        abort();
    }
    Initialise();

    if(ctor==0) {// first file
        ctor = new AcceleratorModelConstructor();
        ctor->NewModel();
    }

    StripHeader((*ifs));

    energy = Pref;
    while((*ifs))
        ReadComponent();
}

AcceleratorModel* MADInterface::GetModel()
{
	assert(ctor);

	if(logFlag && log)
	{
		*log<<endl;
		ctor->ReportStatistics(*log);
		if(inc_sr)
		{
			*log<<"\n\n final energy = "<<energy<<" GeV"<<endl;
		}
	}

	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;
	ctor=0;
	return theModel;
}

AcceleratorModel* MADInterface::ConstructModel ()
{
	if(!ifs)
	{
		MerlinIO::error()<<"MADInterface :: No model file defined!"<<endl;
		abort();
	}

	if(ctor!=0)
	{
		delete ctor;
	}

	ctor = new AcceleratorModelConstructor();
	double z=0.0;
	ctor->NewModel();

	// reset the stream pointer
	// ifs.seekg(0);
	StripHeader((*ifs));

	while(*ifs)
	{
		z+=ReadComponent();
	}

	if(logFlag && log)
	{
		*log<<endl;
		ctor->ReportStatistics(*log);
		*log<<"\nARC distance from MAD file: "<<z<<endl;
		if(inc_sr)
		{
			*log<<"\n\n final energy = "<<energy<<" GeV"<<endl;
		}
	}

	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;
	ctor=0;
	return theModel;
}

void MADInterface::IgnoreZeroLengthType (const string& madType)
{
	zeroLengths.insert(madType);
}

void MADInterface::TreatTypeAsDrift (const std::string& typestr)
{
	driftTypes.insert(typestr);
}

void MADInterface::ConstructNewFrame (const string& name)
{
	SequenceFrame* newFrame;
	if(name[1]!='_')
	{
		if(honMadStructs)
		{
			newFrame = new SequenceFrame(name);
		}
		else
		{
			return;
		}
	}

	else
	{
	switch( name[0] )
	{
		case 'F':
		newFrame = new SequenceFrame(name.substr(2));
		break;

		case 'S':
		newFrame = new SimpleMount(name.substr(2));
		break;

		case 'G':
		//newFrame = new SimpleMount(name.substr(2));
		newFrame = new GirderMount(name.substr(2));
		break;

		case 'M':
		newFrame = new MagnetMover(name.substr(2));
		break;

	        default:
			MERLIN_ERR<<"Unknown frame character: "<<name<<endl;
			abort();
			break;
	}
	}
	ctor->NewFrame(newFrame);
	if(log)
	{
		Log(newFrame->GetName()+" BEGIN",ctor->GetCurrentFrameDepth(),*log);
	}
}//End ConstructNewFrame

void MADInterface::EndFrame (const string& name)
{
	if((!honMadStructs) && name[1]!='_')
	{
		return;
	}

	SequenceFrame& currentFrame = ctor->GetCurrentFrame();

#ifndef NDEBUG
	if(honMadStructs)
	{
		assert(name==currentFrame.GetName());
	}
	else
	{
		assert(name.substr(2)==currentFrame.GetName());
	}
#endif

	if(log)
	{
		Log(currentFrame.GetName()+" END",ctor->GetCurrentFrameDepth(),*log);
	}

	ctor->EndFrame();
}

double MADInterface::ReadComponent ()
{
#define  _READ(value) if(!((*ifs)>>value)) return 0;

	string name,type,aptype;
	double len,ks,angle,e1,e2,k1,k2,k3,h,tilt;
	_READ(name);
	_READ(type);

	prmMap->ReadRow((*ifs));

	name=StripQuotes(name);
	type=StripQuotes(type);

	AcceleratorComponent *component;
	double brho = energy/eV/SpeedOfLight;

	if(prmMap->has_type)
	{
		_READ(aptype);
		aptype=StripQuotes(aptype);
	}
	else if(incApertures)
	{
		MerlinIO::warning()<<"No aperture information. Apertures will not be constructed"<<endl;
		incApertures=false;
	}

	try
	{

        if(driftTypes.find(type)!=driftTypes.end())
	{
		MerlinIO::warning()<<"Treating "<<type<<" as Drift"<<endl;
		type="DRIFT";
        }

        // get the 'standard' parameters
        len = prmMap->GetParameter("L");
        tilt = prmMap->GetParameter("TILT",false);

        if(len==0 && zeroLengths.find(type)!=zeroLengths.end())
        {
		MerlinIO::warning()<<"Ignoring zero length "<<type<<endl;
		return 0;
        }

	if(type=="VKICKER")
	{
		type="YCOR";
	}
	else if(type=="HKICKER")
	{
		type="XCOR";
	}
	else if(type =="KICKER")
	{
		type="DRIFT";
	}
	else if(type =="RFCAVITY")
	{
		//type="DRIFT";
		type="RFCAVITY";
	}
	else if(type=="LCAV")
	{
		type="RFCAVITY";
	}
	else if(type=="RCOLLIMATOR")    // added by Adriana Bungau, 26 October 2006
	{
		type="SPOILER";
	}
	else if(type=="ECOLLIMATOR")    // added by Adriana Bungau, 26 October 2006
	{
		type="SPOILER";
	}

	/******************************************************************
	**
	**	Addition of missing elements in V6.503 LHC "as built" optics
	**
	******************************************************************/

	else if(type=="PLACEHOLDER")
	{
		type="DRIFT";
	}

	else if(type=="TKICKER")//Seems to be standard kickers?
	{
		//type="KICKER"
		type="DRIFT";
	}

	//End of additions - JM

        if(type=="RBEND")
        {
		  if((prmMap->GetParameter("K0L"))!=0.0) type="SBEND";
	}

	if(type=="MULTIPOLE")
	{
		if((prmMap->GetParameter("K0L"))!=0.0) type="SBEND";
		else if((prmMap->GetParameter("K1L"))!=0.0) type="QUADRUPOLE";
		else if((prmMap->GetParameter("K2L"))!=0.0) type="SEXTUPOLE";
		else if((prmMap->GetParameter("K3L"))!=0.0) type="OCTUPOLE";
		else if((prmMap->GetParameter("K4L"))!=0.0) type="DECAPOLE";
		else type="DRIFT";
	}

	if(type=="DRIFT") {
            Drift* aDrift = new Drift(name,len);
            ctor->AppendComponent(*aDrift);
            component=aDrift;
        }
        else if(type=="SPOILER")      // modified by R.Barlow, 30 October 2006
        {
        // double X0 =prmMap->GetParameter("KS"); // cheat! use KS column for radiation length
	Spoiler* aSpoiler = new Spoiler(name,len);
         //|| name[2]=='I' || name[2]=='L'|| name[2]=='P' || name[2]=='S' || name[2]=='D'

	if(name[2]=='T'|| name[2]=='I' || name[2]=='L'|| name[2]=='P' || name[2]=='S' || name[2]=='D')
	{
		aSpoiler->scatter_at_this_spoiler =true;
	}
           //if(apname=="~") apname="TCTH.4L2.B1"

	if(aptype=="~") aptype="T99X99Y99M99";
	double w,h,t;
	int m =0;
	int it=aptype.find("T");
	int ix=aptype.find("X");
	int iy=aptype.find("Y");
	int im=aptype.find("M");
	int ilen=aptype.length();

	//cout << "spoiler: " << aptype << endl;

         // cout<<aptype.substr(it+1,ilen)<<" "<<aptype.substr(ix+1,iy-1)<<" "<<aptype.substr(iy+1,ilen)<<endl;
	istringstream istr(aptype.substr(ix+1,iy-1));
	istr >> w;

	/*
	material type:
	00 Be
	01 C
	02 Al
	03 Cu
	04 W
	05 Pb
	*/
	if (im > 0)
	{
		istringstream mstr(aptype.substr(im+1,ilen));
		mstr>>m;
	}
	else
	{
		im = ilen+1; // stop bad things happening if no material
	}
	

	double ww = w*0.001; // width in meters (AB - Sept 2007)
	istringstream jstr(aptype.substr(iy+1,im-1));
	jstr>>h;
	double hh = h*0.001;//width in meters (AB - Sept 2007)
	istringstream kstr(aptype.substr(it+1,ix-1));
	kstr>>t;
	double tt = t;// tilt in rad

	TiltedAperture* app=new TiltedAperture(ww,hh,tt,m);
	aSpoiler->SetAperture(app);
	ctor->AppendComponent(*aSpoiler);
	component=aSpoiler;
	double X0 = app->X0;
	double conductivity = app->sigma;
	double aperture_size = ww;

	//Collimation only will take place on one axis
	if (hh < ww){aperture_size = hh;} //set to smallest out of hh or ww
	ResistivePotential* resWake =  new ResistivePotential(1,conductivity,0.5*aperture_size,len*meter,"Data/table");
	aSpoiler->SetWakePotentials(resWake);
        } //End Spoiler

	//Magnets

        else if(type=="QUADRUPOLE")
        {
		k1=prmMap->GetParameter("K1L");
		Quadrupole* quad = new Quadrupole(name,len,brho*k1/len);
		ctor->AppendComponent(*quad);
		component=quad;
	}
        else if(type=="SKEWQUAD")
        {
		k1=prmMap->GetParameter("K1L");
		SkewQuadrupole* quad = new SkewQuadrupole(name,len,brho*k1/len);
		ctor->AppendComponent(*quad);
		component=quad;
        }
	else if(type=="SOLENOID")
	{
		ks=prmMap->GetParameter("KS");
		Solenoid* aSolenoid = new Solenoid(name,len,brho*ks/len);
		ctor->AppendComponent(*aSolenoid);
		component=aSolenoid;
	}


	//Dipole Bending Magnets
        else if(type=="SBEND")
        {
		angle=prmMap->GetParameter("K0L");
		k1   =prmMap->GetParameter("K1L");
		h = angle/len;
		SectorBend* bend = new SectorBend(name,len,h,brho*h);

		if(k1!=0)  // mixed function dipole
		{
			bend->SetB1(brho*k1/len);
		}

		e1   =prmMap->GetParameter("E1");
		e2   =prmMap->GetParameter("E2");

		if(e1!=0 || e2!=0)
		{
			if(e1==e2)
			{
				bend->SetPoleFaceInfo(new SectorBend::PoleFace(e1));
			}

			else
			{
				SectorBend::PoleFace* pf1 = e1!=0 ? new SectorBend::PoleFace(e1) : 0;
				SectorBend::PoleFace* pf2 = e2!=0 ? new SectorBend::PoleFace(e2) : 0;
				bend->SetPoleFaceInfo(pf1,pf2);
			}
		}

		if(tilt!=0)
		(*bend).GetGeometry().SetTilt(tilt);

		// check for synchrotron radiation
		if(inc_sr)
		{
			energy -= SRdE(h,len,energy);
		}

		ctor->AppendComponent(*bend);
		component=bend;
	} //End SBEND

	else if(type=="SEXTUPOLE")
	{
		k2=prmMap->GetParameter("K2L");
		Sextupole* sx = new Sextupole(name,len,brho*k2/len);
		ctor->AppendComponent(*sx);
		component=sx;
        }
        else if(type=="OCTUPOLE")
        {
		k3=prmMap->GetParameter("K3L");
		Octupole* oct = new Octupole(name,len,brho*k3/len);
		ctor->AppendComponent(*oct);
		component=oct;
        }
	else if(type=="SKEWSEXT")
	{
		k2=prmMap->GetParameter("K2L");
		SkewSextupole* sx = new SkewSextupole(name,len,brho*k2/len);
		ctor->AppendComponent(*sx);
		component=sx;
        }
        else if(type=="YCOR")
	{
		YCor* yc = new YCor(name,len);
		ctor->AppendComponent(*yc);
		component=yc;
        }
	else if(type=="XCOR")
	{
		XCor* xc = new XCor(name,len);
		ctor->AppendComponent(*xc);
		component=xc;
	}

        else if(type=="RFCAVITY")
        {
		cout << "Found RF cavity\t";
		// Here we assume an SW cavity
		double freq=prmMap->GetParameter("FREQ");
		double phase=prmMap->GetParameter("LAG");
		double volts=prmMap->GetParameter("VOLT");

		// standing wave cavities need an exact integer of half-wavelengths
		freq*=MHz;
		double lambdaOver2 = SpeedOfLight/freq/2;
		int ncells = Round(len/lambdaOver2);
		double len1 = ncells*lambdaOver2;

		cout << "f: " << freq << "\tV: " << volts << "\tncells: " << ncells << "\tWavelength/2: " << lambdaOver2 << "\tLength: " << len1 << endl;

		// adjust phase for cosine-like field
		phase = twoPi*(phase-0.25);

		if((len1/len-1) > 0.001)
		{
			MerlinIO::error() << "SW cavity length not valid ";
			MerlinIO::error() << '('<<len<<", "<<len1<<')'<<endl;
		}

		SWRFStructure* rfsctruct = new SWRFStructure(name,ncells,freq,volts*MV/len,phase);
		ctor->AppendComponent(*rfsctruct);
		component=rfsctruct;
	}//End RFCAVITY

        else if(type=="CRABRF")
	{
		TransverseRFStructure* crabcav = new TransverseRFStructure(name,len,0,0);
		ctor->AppendComponent(*crabcav);
		component=crabcav;
        }

        else if(type=="MONITOR")
        {
		if(name.substr(0,4)=="BPM_")
		{
			BPM* bpm = new BPM("BPM"+name.substr(4));
			ctor->AppendComponent(*bpm);
			component=bpm;
		}
		else if(name.substr(0,3)=="WS_")
		{
			RMSProfileMonitor* ws = new RMSProfileMonitor("WS"+name.substr(3));
			ctor->AppendComponent(*ws);
			component=ws;
		}
		else
		{
			MERLIN_WARN<<"unknown monitor type: "<<name<<" defaulting to BPM"<<endl;
			BPM* bpm = new BPM(name);
			ctor->AppendComponent(*bpm);
			component=bpm;
		}
	}
        else if(type=="LINE")
        {
		if(!flatLattice)
		{
			if(!frameStack.empty() && name == frameStack.top())
			{
				frameStack.pop();
				EndFrame(name);
			}

			else
			{
				frameStack.push(name);
				ConstructNewFrame(name);
			}
		}
		component=0;
	}

	else if(type=="MATRIX") // just ignore for now.
	{
		component=0;
	}

	else if(type=="SROT")
	{
		ctor->AppendComponentFrame(ConstructSrot(prmMap->GetParameter("L"),name));
		component=0;
	}
        else
	{
		MERLIN_ERR<<"ERROR: undefined type: "<<type<<endl;
		abort();
	}

	if(component && log)
	{
		Log(component->GetQualifiedName(),ctor->GetCurrentFrameDepth(),*log);
	}

	if(component && incApertures)
	{
		component->SetAperture(ConstructAperture(aptype));
	}

	}//End of try block

	catch(MADKeyMap::bad_key)
	{
		MerlinIO::error()<<"L not present in table"<<endl;
		abort();
	}

	return component ? component->GetLength() : 0.0;
}
