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

#include "IO/MerlinIO.h"
#include "NumericalUtils/utils.h"

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/Apertures/CollimatorAperture.h"
#include "AcceleratorModel/Apertures/RectEllipseAperture.h"
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"
#include "AcceleratorModel/Frames/SequenceFrame.h"
#include "AcceleratorModel/Supports/SupportStructure.h"
#include "AcceleratorModel/Supports/MagnetMover.h"

#include "Collimators/ResistiveWakePotentials.h"

#include "NumericalUtils/PhysicalConstants.h"

#include "MADInterface/ConstructSrot.h"
#include "MADInterface/MADKeyMap.h"
#include "MADInterface/MADInterface.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

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
	char buffer[1024];
	char c;
	while(true && is)
	{
		is.get(c);
		if(c=='*'||c=='$'||c=='@')
		{
			is.getline(buffer,1024);
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

//Aperture* ConstructAperture(const string& apstr)
Aperture* ConstructAperture(const double& ap_type, MADKeyMap* prmMap)
{
	Aperture* ap;
//	double w,h,d,t,r,a,b;
//	string::size_type n,n2;
	double w,h,r,a,b;

	/*
	Aperture types - from MADX: http://mad.web.cern.ch/mad/Introduction/aperture.html
	CIRCLE		1
	ELLIPSE		2
	RECTANGLE	3
	LHCSCREEN	4
	MARGUERITE	5
	RECTELLIPSE	6
	RACETRACK	7
	NONE		0

	aper_# means for all apertypes but racetrack:
	aper_1 = half width rectangle
	aper_2 = half height rectangle
	aper_3 = half horizontal axis ellipse (or radius if circle)
	aper_4 = half vertical axis ellipse

	For racetrack, the aperture parameters will have the same meaning as the tolerances:
	aper_1 and xtol = horizontal displacement of radial part
	aper_2 and ytol = vertical displacement of radial part
	aper_3 and rtol = radius
	aper_4 = not used 
	*/


	//Circle
	if(ap_type == 1)
	{
		r = prmMap->GetParameter("APER_3");
		if(r == 0.0)
		{
			//Zero radius, disable the aperture
			ap = 0;
		}
		else
		{
			//We have a non-zero radius, create the aperture
			ap = new CircularAperture(r);
		}
	}

	//RECTANGLE
	else if(ap_type == 3)
	{
		w = prmMap->GetParameter("APER_1");	//half width rectangle
		h = prmMap->GetParameter("APER_2");	//half height rectangle

		if (w == 0.0 || h == 0.0)
		{
			ap = 0;
		}
		else
		{
			ap = new RectangularAperture(w,h);
		}
	}

	//RECTELLIPSE
	else if(ap_type == 6)
	{
		//FIXME
		w = prmMap->GetParameter("APER_1");	//half width rectangle
		h = prmMap->GetParameter("APER_2");	//half height rectangle
		a = prmMap->GetParameter("APER_3");	//half horizontal axis ellipse
		b = prmMap->GetParameter("APER_4");	//half vertical axis ellipse

		if (w == 0.0 || h == 0.0 || a == 0.0 || b == 0.0)
		{
			ap = 0;
		}
		else
		{
			ap = new RectEllipseAperture (w, h, a, b);
			//cout << endl << endl << w <<"\t" << h << endl;
		}
	}

	//NONE
	else if(ap_type == 0.0)
	{
		//No aperture defined, disable aperture for this component
		ap = 0;
	}

	else
	{
		MERLIN_ERR << "WARNING: unknown aperture definition ("<<ap_type<<") ignored"<<endl;
		ap=0;
	}

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

} // Namespace end

// Class MADInterface
/*
MADInterface::MADInterface (const std::string& madFileName, double P0)
        : energy(P0),ifs(madFileName.empty() ? 0 : new ifstream(madFileName.c_str())),
        log(MerlinIO::std_out),logFlag(false),flatLattice(false),honMadStructs(false),
        incApertures(true),inc_sr(false),ctor(0),prmMap(0),collimator_db(NULL),z(0)
*/
MADInterface::MADInterface (const std::string& madFileName, double P0)
        : energy(P0),ifs(madFileName.empty() ? 0 : new ifstream(madFileName.c_str())),
        log(MerlinIO::std_out),logFlag(false),flatLattice(false),honMadStructs(false),
        incApertures(true),inc_sr(false),ctor(0),prmMap(0),z(0)
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
//	TreatTypeAsDrift("MARKER");	// merlin bug!
	TreatTypeAsDrift("INSTRUMENT"); // merlin bug!

	//Addition of missing elements in V6.503 LHC "as built" optics
	TreatTypeAsDrift("PLACEHOLDER"); // placeholders for extra upgrade components etc (LHC)	
	
	TreatTypeAsDrift("TKICKER");	// merlin bug! - transverse dampers, injection + extraction kickers + friends.
	//TreatTypeAsDrift("RFCAVITY");	// merlin bug! - Fix tracking with zero cavity voltage.

	IgnoreZeroLengthType("RCOLLIMATOR");
}

MADInterface::~MADInterface(){
	if(ctor)
		delete ctor;
	if(prmMap)
		delete prmMap;
	if(ifs)
		ifs->close();
		delete ifs;
}

void MADInterface::Initialise()
{
	if(prmMap!=0)
	{
		delete prmMap;
	}

	string s;
	char c;
	bool tfs = 0;
	//This just grabs the top line, we really want the line that starts with *, aka the one with the headers.
	//The following checks were for a file that is not valid TFS format, changing to check for what MADX actually produces - JM
	//getline((*ifs),s);
	//check_column_heading((*ifs),"*");
	//check_column_heading((*ifs),"NAME");
	//check_column_heading((*ifs),"KEYWORD");
	//Is MAD8 different?

	//check_column_heading((*ifs),"@");
	//check_column_heading((*ifs),"NAME");	
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
	{
		delete ifs;
	}

	ifs =  new ifstream(fname.c_str());

	if(!(*ifs))
	{
		MERLIN_ERR << "ERROR opening file " << fname << endl;
		delete ifs;
		abort();
	}

	Initialise();

	if(ctor==0)
	{
		// first file
		ctor = new AcceleratorModelConstructor();
		ctor->NewModel();
	}

	StripHeader((*ifs));

	energy = Pref;
	while((*ifs))
	{
		ReadComponent();
	}
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
			*log << endl << endl << "final energy = " << energy << " GeV" << endl;
		}
	}

	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;
	ctor=0;
	return theModel;
}

AcceleratorModel* MADInterface::ConstructModel()
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
	//double z=0.0;
	ctor->NewModel();

	// reset the stream pointer
	// ifs.seekg(0);
	StripHeader((*ifs));


	//cout << "Name\tType\tS" << endl;

	//Main component read in loop
	while((*ifs).good())
	{
		z+=ReadComponent();
	}

	if(logFlag && log)
	{
		*log << endl;
		ctor->ReportStatistics(*log);
		*log << endl << "ARC distance from MAD file: " << z << endl;
		if(inc_sr)
		{
			*log << endl << endl << "final energy = " << energy << " GeV" << endl;
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
				MERLIN_ERR << "Unknown frame character: " << name << endl;
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

	string name,type,aptype,parent,aperture;
	double len,ks,angle,e1,e2,k1,k2,k3,h,tilt;
	_READ(name);
	_READ(type);

	prmMap->ReadRow((*ifs));

	name=StripQuotes(name);
	type=StripQuotes(type);

	//cout << name << "\t" << type << "\t" << prmMap->GetParameter("S") << endl;

	AcceleratorComponent *component = NULL;
	double brho = energy/eV/SpeedOfLight;

	if(prmMap->has_type)
	{
		_READ(aptype);
		aptype=StripQuotes(aptype);
	}

	//Do we want to build apertures, and do we have the required information required?
	if(incApertures && !prmMap->has_apertype)
	{
		//Could not find aperture information and building of apertures was requested, will disable building of apertures.
		MerlinIO::warning() << "No aperture information. Apertures will not be constructed" << endl;
		incApertures=false;
	}

	try
	{

	if(driftTypes.find(type)!=driftTypes.end())
	{
		#ifndef NDEBUG
		MerlinIO::warning() << "Treating " << type << " as Drift" << endl;
		#endif

		type="DRIFT";
	}

        // get the 'standard' parameters
        len = prmMap->GetParameter("L");
        tilt = prmMap->GetParameter("TILT",false);

        if(len==0 && zeroLengths.find(type)!=zeroLengths.end())
        {
		MerlinIO::warning() << "Ignoring zero length " << type << ": " << name << endl;
		return 0;
        }
	else if(type =="KICKER")
	{
		type="DRIFT";
	}
	else if(type =="TKICKER")
	{
		type="DRIFT";
	}
/*	else if(type=="VKICKER")
	{
		if (prmMap->GetParameter("VKICK") == 0)
		{
			type="DRIFT";
		}
	}
	else if(type=="HKICKER")
	{
		if (prmMap->GetParameter("HKICK") == 0)
		{
			type="DRIFT";
		}
	}
*/
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
		type="COLLIMATOR";
	}
	else if(type=="ECOLLIMATOR")    // added by Adriana Bungau, 26 October 2006
	{
		type="COLLIMATOR";
	}

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

	if(type=="DRIFT")
	{
		if(len !=0)
		{
		Drift* aDrift = new Drift(name,len);
		ctor->AppendComponent(*aDrift);
		component=aDrift;
		}
		else
		{
			component = 0;
		}
        }
	else if(type=="VKICKER")
	{
		double scale;
		if(len > 0)
		{
			scale = brho/len;
		}
		else	//treat as an integrated length
		{
			scale = brho;
		}

		double kick = prmMap->GetParameter("VKICK");
		//X,Y,tilt
	//	cout << "VKICKER " << name << "\t" << len << "\t" << kick << "\t" << tilt << endl;
		YCor* aKicker = new YCor(name,len,scale*kick);
//		if(tilt!=0)
//		(*aKicker).GetGeometry().SetTilt(tilt);
		ctor->AppendComponent(*aKicker);
		component=aKicker;
	}
	else if(type=="HKICKER")
	{
		double scale;
		if(len > 0)
		{
			scale = brho/len;
		}
		else	//treat as an integrated length
		{
			scale = brho;
		}
		double kick = prmMap->GetParameter("HKICK");
		//X,Y,tilt
	//	cout << "HKICKER " << name << "\t" << len << "\t" << kick << "\t" << tilt << endl;
		XCor* aKicker = new XCor(name,len,-scale*kick);
//		if(tilt!=0)
//		(*aKicker).GetGeometry().SetTilt(tilt);
		ctor->AppendComponent(*aKicker);
		component=aKicker;
	}
	else if(type=="COLLIMATOR")
	{
		/*
		//Check we have collimator info
		if(collimator_db == NULL)
		{
			std::cerr << "Collimator settings not found. Exiting." << std::endl;
			exit(EXIT_FAILURE);
		}

		//Also check that we are not using a zero length collimator
		if(len == 0.0)
		{
			cout << "Rejecting " << name << " due to zero length" << endl;
			return 0;
		}

		//Lets set up some variables for collimator settings
		double collimator_aperture_width;
		double collimator_aperture_height;
		double collimator_aperture_tilt;
		material* collimator_material;
		
		//We now need to find the collimator configuration
		bool have_collimator = false;

		for(unsigned int i=1; i < collimator_db->number_collimators; i++)
		{
			//Time to search for the collimator we are currently using
			if(collimator_db->Collimator[i].name == name)
			{
				#ifndef NDEBUG
				cout << "Found the collimator: " << collimator_db->Collimator[i].name << endl;
				#endif
				have_collimator = true;

				//This is the collimator we are using. Input file should have half gaps.
				//Factor of 2 turns it into full-gaps which the rest of the code expects for now.

				if(!collimator_db->use_sigma)
				{
					collimator_aperture_width = collimator_db->Collimator[i].x_gap*2.0;
					collimator_aperture_height = collimator_db->Collimator[i].y_gap*2.0;
					collimator_aperture_tilt = collimator_db->Collimator[i].tilt;
					collimator_material = collimator_db->Collimator[i].Material;
				}
				else if(collimator_db->use_sigma)
				{
					//Just want to set the distance for lookup after the lattice is completed and then calculate the twiss params
					//The aperture gap will then be calculated.
					collimator_db->Collimator[i].position = z;
				}
				
				//Create a new collimator
				Collimator* aCollimator = new Collimator(name,len);

				//Enable scattering
				aCollimator->scatter_at_this_collimator = true;

				//Will calculate the collimator jaw aperture later; same with any collimator related wakes
				//Use Collimator_Database::Configure_collimators() to do so.
				if(!collimator_db->use_sigma)
				{
					//Create an aperture for the collimator jaws
					CollimatorAperture* app=new CollimatorAperture(collimator_aperture_width,collimator_aperture_height,collimator_aperture_tilt,collimator_material);
		
					//Set the aperture for collimation
					aCollimator->SetAperture(app);
				}

				//Add the component to the accelerator
				ctor->AppendComponent(*aCollimator);
				component=aCollimator;

				//Again, only if we have already specified the aperture
				if(!collimator_db->use_sigma)
				{
					double conductivity = collimator_material->sigma;
					double aperture_size = collimator_aperture_width;

					//Collimation only will take place on one axis
					if (collimator_aperture_height < collimator_aperture_width)
					{
						aperture_size = collimator_aperture_height;
					} //set to smallest out of height or width
				
					//Define the resistive wake for the collimator jaws.
					ResistivePotential* resWake =  new ResistivePotential(1,conductivity,0.5*aperture_size,len*meter,"Data/table");

					//Set the Wake potentials for this collimator
					aCollimator->SetWakePotentials(resWake);
				}
			}

			//We did not find the settings for the collimator in question
		}

		//If we don't have the collimator settings, and hence have left the collimator wide open, then don't enable the wakefields.
		if(have_collimator == false)
		{
			
			//Exit or just default to wide open?
			//Default to wide open.
			#ifndef NDEBUG
			std::cerr << "Collimator " << name << " settings not found in collimator database." << std::endl;
			#endif
			
			//There is no point in wasting cpu time on aperture checks/wakefields with a wide open collimator
			//Lets just treat the element as a drift instead.
			Drift* aDrift = new Drift(name,len);
			ctor->AppendComponent(*aDrift);
			component=aDrift;
		}
		*/


		Collimator* aCollimator = new Collimator(name,len);

		//Add the component to the accelerator
		ctor->AppendComponent(*aCollimator);
		component=aCollimator;
	}//End of Collimators

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
/*			if(e1==e2)
			{
				bend->SetPoleFaceInfo(new SectorBend::PoleFace(e1));
			}

			else
			{
				SectorBend::PoleFace* pf1 = e1!=0 ? new SectorBend::PoleFace(e1) : 0;
				SectorBend::PoleFace* pf2 = e2!=0 ? new SectorBend::PoleFace(e2) : 0;
				bend->SetPoleFaceInfo(pf1,pf2);
			}
*/
			SectorBend::PoleFace* pf1 = e1!=0 ? new SectorBend::PoleFace(e1) : 0;
			SectorBend::PoleFace* pf2 = e2!=0 ? new SectorBend::PoleFace(e2) : 0;
			bend->SetPoleFaceInfo(pf1,pf2);
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
	
	else if(type=="HEL")
	{
		HollowElectronLens* hel = new HollowElectronLens(name, len);
		ctor->AppendComponent(*hel);
		component=hel;
	}

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
			BPM* bpm = new BPM("BPM"+name.substr(4),len);
			ctor->AppendComponent(*bpm);
			component=bpm;
		}
		else if(name.substr(0,3)=="WS_")
		{
			RMSProfileMonitor* ws = new RMSProfileMonitor("WS"+name.substr(3),len);
			ctor->AppendComponent(*ws);
			component=ws;
		}
		else
		{
			#ifndef NDEBUG
			MERLIN_WARN << "Unknown monitor type: "<<name<<" defaulting to BPM" << endl;
			#endif
			BPM* bpm = new BPM(name,len);
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
	else if(type=="MARKER")
	{
		Marker* aMarker = new Marker(name);
		ctor->AppendComponent(*aMarker);
		component=aMarker;
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

	if(component && incApertures && type!="COLLIMATOR")
	{
		component->SetAperture(ConstructAperture(prmMap->GetParameter("APERTYPE"),prmMap));
	}

	}//End of try block

	catch(MADKeyMap::bad_key)
	{
		MerlinIO::error()<<"L not present in table"<<endl;
		abort();
	}
	if(component)
	{
		component->SetComponentLatticePosition(z);
	}

	return component ? component->GetLength() : 0.0;
}
