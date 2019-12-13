/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include "Components.h"
#include "MerlinIO.h"
#include "AcceleratorModelConstructor.h"
#include "SequenceFrame.h"
#include "SupportStructure.h"
#include "MagnetMover.h"
#include "ResistiveWakePotentials.h"
#include "PhysicalConstants.h"
#include "DataTableTFS.h"
#include "ConstructSrot.h"
#include "MADInterface.h"

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

namespace
{

//stack<string> frameStack;

void Log(const string& tag, int depth, ostream& os)
{
	static const char* tab = "----|";
	while(depth--)
	{
		os << tab;
	}
	os << ' ' << tag << endl;
}

}

// Class MADInterface
MADInterface::MADInterface(const string& madFileName, double P0) :
	momentum(P0), filename(madFileName)
{
	infile = make_unique<ifstream>(madFileName);
	ifs = infile.get();
	init();
}

MADInterface::MADInterface(std::istream *in, double P0) :
	momentum(P0), filename("std::istream"), ifs(in)
{
	init();
}

void MADInterface::init()
{
	if(!ifs || !ifs->good())
	{
		MERLIN_ERR << "MADInterface: ERROR opening or reading file " << filename << std::endl;
		throw MerlinException(string("ERROR opening file ") + string(filename));
	}

	log = MerlinIO::std_out;

	//By default, we currently treat the following MAD types as drifts
	TreatTypeAsDrift("INSTRUMENT");
	TreatTypeAsDrift("PLACEHOLDER");
	TreatTypeAsDrift("VMONITOR");
	TreatTypeAsDrift("HMONITOR");
	TreatTypeAsDrift("KICKER");
	TreatTypeAsDrift("TKICKER"); // merlin bug! - transverse dampers, injection + extraction kickers + friends.
	TreatTypeAsDrift("MATRIX");

	IgnoreZeroLengthType("RCOLLIMATOR");
	//IgnoreZeroLengthType("ECOLLIMATOR");
}

MADInterface::~MADInterface()
{
	if(modelconstr)
	{
		delete modelconstr;
	}
}

inline double SRdE(double h, double len, double E)
{
	// Calculation assumes electron, and uses energy rather than momentum
	static const double Cg = 8.85e-05 / twoPi;
	return Cg * pow(E, 4) * h * h * len;
}

AcceleratorModel* MADInterface::ConstructModel()
{
	unique_ptr<DataTable> MADinput;
	try
	{
		MADinput = make_unique<DataTable>(DataTableReaderTFS(ifs).Read());
	}
	catch(const BadFormatException &e)
	{
		MERLIN_ERR << "MADInterface: Error reading " << filename << endl;
		throw e;
	}

	if(modelconstr != nullptr && appendFlag == false)
	{
		delete modelconstr;
	}
	if(modelconstr == nullptr && appendFlag == false)
	{
		modelconstr = new AcceleratorModelConstructor();
	}

	TypeFactory* factory = new TypeFactory();
	double brho = momentum / eV / SpeedOfLight;

	//Loop over all components
	for(auto &MADinputrow : *MADinput)
	{
		string type = MADinputrow.Get_s("KEYWORD");
		double length = MADinputrow.Get_d("L");

		if(length == 0 && zeroLengths.find(type) != zeroLengths.end())
		{
			MerlinIO::warning() << "Ignoring zero length " << type << ": " << MADinputrow.Get_s("NAME") << endl;
			continue;
		}
		TypeOverrides(MADinputrow);

		if(type == "LINE")
		{
			if(!flatLattice)
			{
				const string& name = MADinputrow.Get_s("NAME");
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
			continue;
		}
		else if(type == "SROT")
		{
			modelconstr->AppendComponentFrame(ConstructSrot(length, MADinputrow.Get_s("NAME")));
			continue;
		}

		//Determine multipole type by parameters
		vector<AcceleratorComponent*> components = factory->GetInstance(MADinputrow, brho);

		if(inc_sr && (type == "SBEND" || type == "RBEND"))
		{
			momentum -= SRdE(MADinputrow.Get_d("ANGLE") / length, length, momentum);
			brho = momentum / eV / SpeedOfLight;
		}

		for(auto component : components)
		{
			modelconstr->AppendComponent(*component);
			component->SetComponentLatticePosition(z);
			z += component->GetLength();
		}
	} //End for loop

	if(logFlag && log)
	{
		*log << endl;
		modelconstr->ReportStatistics(*log);
		*log << endl << "ARC distance from MAD file: " << z << endl;

		if(inc_sr)
		{
			*log << endl << endl << "final momentum = " << momentum << " GeV" << endl;
		}
	}

	AcceleratorModel* theModel = modelconstr->GetModel();
	delete modelconstr;
	delete factory;
	modelconstr = nullptr;
	return theModel;
}

void MADInterface::TypeOverrides(DataTableRow& MADinputrow)
{
	string keyword = MADinputrow.Get_s("KEYWORD");
	if(driftTypes.find(keyword) != driftTypes.end())
		MADinputrow.Set_s("KEYWORD", "DRIFT");
	if(keyword == "LCAV")
		MADinputrow.Set_s("KEYWORD", "RFCAVITY");
	if(keyword == "RCOLLIMATOR" || keyword == "ECOLLIMATOR")
		MADinputrow.Set_s("KEYWORD", "COLLIMATOR");
	if(keyword == "RBEND" && MADinputrow.Get_d("K0L"))
		MADinputrow.Set_s("KEYWORD", "SBEND");
	if(single_cell_rf && MADinputrow.Get_s("KEYWORD") == "RFCAVITY")
		MADinputrow.Set_s("KEYWORD", "RFCAVITY_SingleCell");
}

string MADInterface::GetMutipoleType(DataTableRow& MADinputrow)
{
	if(!MADinputrow.Get_d("K0L"))
		return "SBEND";
	if(!MADinputrow.Get_d("K1L"))
		return "QUADRUPOLE";
	if(!MADinputrow.Get_d("K2L"))
		return "SEXTUPOLE";
	if(!MADinputrow.Get_d("K3L"))
		return "OCTUPOLE";
	if(!MADinputrow.Get_d("K4L"))
		return "DECAPOLE";
	return "DRIFT";
}

void MADInterface::ConstructNewFrame(const string& name)
{
	SequenceFrame* newFrame;
	if(name[1] != '_')
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
		switch(name[0])
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
	modelconstr->NewFrame(newFrame);

	if(log)
	{
		Log(newFrame->GetName() + " BEGIN", modelconstr->GetCurrentFrameDepth(), *log);
	}
} //End ConstructNewFrame

void MADInterface::EndFrame(const string& name)
{
	if((!honMadStructs) && name[1] != '_')
	{
		return;
	}

	SequenceFrame& currentFrame = modelconstr->GetCurrentFrame();

	modelconstr->EndFrame();

	if(log)
	{
		Log(currentFrame.GetName() + " END", modelconstr->GetCurrentFrameDepth(), *log);
	}
}

void MADInterface::AppendModel(const string& fname, double Pref)
{
	filename = fname;

	infile = make_unique<ifstream>(fname);
	ifs = infile.get();

	if(!(*ifs))
	{
		MERLIN_ERR << "ERROR opening file " << fname << endl;
		abort();
	}

	if(modelconstr == nullptr)
	{
		// first file
		modelconstr = new AcceleratorModelConstructor();
	}

	momentum = Pref;

	appendFlag = 1;
	ConstructModel();
}

AcceleratorModel* MADInterface::GetModel()
{
	assert(modelconstr);

	if(logFlag && log)
	{
		*log << endl;
		modelconstr->ReportStatistics(*log);
		if(inc_sr)
		{
			*log << endl << endl << "final momentum = " << momentum << " GeV" << endl;
		}
	}

	AcceleratorModel* theModel = modelconstr->GetModel();
	delete modelconstr;
	modelconstr = nullptr;
	return theModel;
}

AcceleratorModelConstructor* MADInterface::GetModelConstructor()
{
	assert(modelconstr);
	return modelconstr;
}

void MADInterface::IgnoreZeroLengthType(const string& madType)
{
	zeroLengths.insert(madType);
}

void MADInterface::TreatTypeAsDrift(const std::string& typestr)
{
	driftTypes.insert(typestr);
}

vector<AcceleratorComponent*> DriftComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	if(length != 0)
		return {new Drift(name, length)};
	else
		return {};
}

vector<AcceleratorComponent*> RBendComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double angle = MADinputrow.Get_d("ANGLE");
	double k1l = MADinputrow.Get_d("K1L");
	double tilt = MADinputrow.Get_d("TILT");
	double h = angle / length;

	SectorBend* bend = new SectorBend(name, length, h, brho * h);

	if(k1l)
		bend->SetB1(brho * k1l / length);

	double e1 = MADinputrow.Get_d("E1");
	double e2 = MADinputrow.Get_d("E2");

	if(e1 != 0 || e2 != 0)
	{
		SectorBend::PoleFace* pf1 = e1 != 0 ? new SectorBend::PoleFace(e1) : nullptr;
		SectorBend::PoleFace* pf2 = e2 != 0 ? new SectorBend::PoleFace(e2) : nullptr;
		bend->SetPoleFaceInfo(pf1, pf2);
	}
	bend->GetGeometry().SetTilt(tilt);

	return {bend};
}

vector<AcceleratorComponent*> SBendComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double angle = MADinputrow.Get_d("ANGLE");
	double k1l = MADinputrow.Get_d("K1L");
	double tilt = MADinputrow.Get_d("TILT");
	double h = angle / length;

	SectorBend* bend = new SectorBend(name, length, h, brho * h);

	if(k1l)
		bend->SetB1(brho * k1l / length);

	double e1 = MADinputrow.Get_d("E1");
	double e2 = MADinputrow.Get_d("E2");

	if(e1 || e2)
	{
		SectorBend::PoleFace* pf1 = e1 != 0 ? new SectorBend::PoleFace(e1) : nullptr;
		SectorBend::PoleFace* pf2 = e2 != 0 ? new SectorBend::PoleFace(e2) : nullptr;
		bend->SetPoleFaceInfo(pf1, pf2);
	}
	if(tilt)
		bend->GetGeometry().SetTilt(tilt);

	return {bend};
}

vector<AcceleratorComponent*> QuadrupoleComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double k1l = MADinputrow.Get_d("K1L");

	return {new Quadrupole(name, length, brho * k1l / length)};
}

vector<AcceleratorComponent*> SkewQuadrupoleComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double k1l = MADinputrow.Get_d("K1L");

	return {new SkewQuadrupole(name, length, brho * k1l / length)};
}

vector<AcceleratorComponent*> SextupoleComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double k2l = MADinputrow.Get_d("K2L");

	return {new Sextupole(name, length, brho * k2l / length)};
}

vector<AcceleratorComponent*> SkewSextupoleComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double k2l = MADinputrow.Get_d("K2L");

	return {new SkewSextupole(name, length, brho * k2l / length)};
}

vector<AcceleratorComponent*> OctupoleComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double k3l = MADinputrow.Get_d("K3L");

	return {new Octupole(name, length, brho * k3l / length)};
}

vector<AcceleratorComponent*> YCorComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new YCor(name, length)};
}

vector<AcceleratorComponent*> XCorComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new XCor(name, length)};
}

vector<AcceleratorComponent*> VKickerComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double kick = MADinputrow.Get_d("VKICK");
	double scale;
	if(length > 0)
		scale = brho / length;
	else
		scale = brho;

	return {new YCor(name, length, scale * kick)};
}

vector<AcceleratorComponent*> HKickerComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double kick = MADinputrow.Get_d("HKICK");
	double scale;
	if(length > 0)
		scale = brho / length;
	else
		scale = brho;

	return {new XCor(name, length, -scale * kick)};
}

vector<AcceleratorComponent*> SolenoidComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	double ks = MADinputrow.Get_d("KS");

	return {new Solenoid(name, length, brho * ks / length)};
}

vector<AcceleratorComponent*> RFCavityComponentSingleCell::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	// Here we assume an SW cavity
	double freq = MADinputrow.Get_d("FREQ");
	double phase = MADinputrow.Get_d("LAG");
	double volts = MADinputrow.Get_d("VOLT");
	// standing wave cavities need an exact integer of half-wavelengths
	freq *= MHz;
	double lambdaOver2 = SpeedOfLight / freq / 2;
	int ncells = round(length / lambdaOver2);
	double len1 = ncells * lambdaOver2;

	// adjust phase for cosine-like field
	phase = twoPi * (phase - 0.25);

	double rfcav_len = lambdaOver2;
	double drift_len = length - lambdaOver2;
	ncells = round(rfcav_len / lambdaOver2);

	if((rfcav_len / length - 1) > 0.001)
	{
		MerlinIO::error() << "SW cavity length not valid ";
		MerlinIO::error() << '(' << length << ", " << len1 << ')' << endl;
	}

	SWRFStructure* rfstruct = new SWRFStructure(name, ncells, freq, volts * MV / rfcav_len, phase);

	string drift_name = "Drift_" + name;
	Drift* rf_drift = new Drift(drift_name, drift_len);

	return {rfstruct, rf_drift};

}

vector<AcceleratorComponent*> RFCavityComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");
	// Here we assume an SW cavity
	double freq = MADinputrow.Get_d("FREQ");
	double phase = MADinputrow.Get_d("LAG");
	double volts = MADinputrow.Get_d("VOLT");
	// standing wave cavities need an exact integer of half-wavelengths
	freq *= MHz;
	double lambdaOver2 = SpeedOfLight / freq / 2;
	int ncells = round(length / lambdaOver2);
	double len1 = ncells * lambdaOver2;

	// adjust phase for cosine-like field
	phase = twoPi * (phase - 0.25);

	if(((len1 / length) - 1) > 0.001)
	{
		MerlinIO::error() << "SW cavity length not valid ";
		MerlinIO::error() << '(' << length << ", " << len1 << ')' << endl;
	}

	return {new SWRFStructure(name, ncells, freq, volts * MV / length, phase)};

}

vector<AcceleratorComponent*> CrabMarkerComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	double mux = MADinputrow.Get_d("MUX");
	double muy = MADinputrow.Get_d("MUY");
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new CrabMarker(name, length, mux, muy)};
}

vector<AcceleratorComponent*> CrabRFComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new TransverseRFStructure(name, length, 0, 0)};
}

vector<AcceleratorComponent*> CollimatorComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new Collimator(name, length)};
}

vector<AcceleratorComponent*> HELComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	return {new HollowElectronLens(name, length, 0, 0, 0, 0, 0)};
}

vector<AcceleratorComponent*> MonitorComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");
	double length = MADinputrow.Get_d("L");

	if(name.substr(0, 2) == "WS")
		return {new RMSProfileMonitor(name, length)};
	else
		return {new BPM(name, length)};
}

vector<AcceleratorComponent*> MarkerComponent::GetInstance(DataTableRow& MADinputrow, double brho)
{
	const string& name = MADinputrow.Get_s("NAME");

	return {new Marker(name)};
}

vector<AcceleratorComponent*> TypeFactory::GetInstance(DataTableRow& MADinputrow, double brho)
{
	string type = MADinputrow.Get_s("KEYWORD");
	map<string, getTypeFunc>::iterator itr = componentTypes.find(type);
	if(itr != componentTypes.end())
	{
		return (*itr->second)(MADinputrow, brho);
	}
	return {};
}

TypeFactoryInit::TypeFactoryInit()
{
	TypeFactory::componentTypes["DRIFT"] = &DriftComponent::GetInstance;
	TypeFactory::componentTypes["RBEND"] = &RBendComponent::GetInstance;
	TypeFactory::componentTypes["SBEND"] = &SBendComponent::GetInstance;
	TypeFactory::componentTypes["QUADRUPOLE"] = &QuadrupoleComponent::GetInstance;
	TypeFactory::componentTypes["SKEWQUAD"] = &SkewQuadrupoleComponent::GetInstance;
	TypeFactory::componentTypes["SEXTUPOLE"] = &SextupoleComponent::GetInstance;
	TypeFactory::componentTypes["SKEWSEXT"] = &SkewSextupoleComponent::GetInstance;
	TypeFactory::componentTypes["OCTUPOLE"] = &OctupoleComponent::GetInstance;
	TypeFactory::componentTypes["YCOR"] = &YCorComponent::GetInstance;
	TypeFactory::componentTypes["XCOR"] = &XCorComponent::GetInstance;
	TypeFactory::componentTypes["VKICKER"] = &VKickerComponent::GetInstance;
	TypeFactory::componentTypes["HKICKER"] = &HKickerComponent::GetInstance;
	TypeFactory::componentTypes["SOLENOID"] = &SolenoidComponent::GetInstance;
	TypeFactory::componentTypes["RFCAVITY"] = &RFCavityComponent::GetInstance;
	TypeFactory::componentTypes["RFCAVITY_SingleCell"] = &RFCavityComponentSingleCell::GetInstance;
	TypeFactory::componentTypes["CRABMARKER"] = &CrabMarkerComponent::GetInstance;
	TypeFactory::componentTypes["CRABRF"] = &CrabRFComponent::GetInstance;
	TypeFactory::componentTypes["COLLIMATOR"] = &CollimatorComponent::GetInstance;
	TypeFactory::componentTypes["HEL"] = &HELComponent::GetInstance;
	TypeFactory::componentTypes["MONITOR"] = &MonitorComponent::GetInstance;
	TypeFactory::componentTypes["MARKER"] = &MarkerComponent::GetInstance;
}

map<string, getTypeFunc> TypeFactory::componentTypes;
TypeFactoryInit TypeFactoryInit::init;
