/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include "Components.h"
#include "MerlinIO.h"
#include "utils.h"
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
	energy(P0), inc_sr(false), flatLattice(false),  z(0), single_cell_rf(false), filename(madFileName), ifs(
		madFileName.empty() ?
		nullptr : new ifstream(madFileName.c_str())), log(
		MerlinIO::std_out), logFlag(false), honMadStructs(false), appendFlag(false), modelconstr(
		nullptr)
{
	if(ifs)
	{
		if(!(*ifs))
		{
			MERLIN_ERR << "ERROR opening file " << madFileName << std::endl;
			throw MerlinException(string("ERROR opening file ") + string(madFileName));
		}
	}
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
	if(ifs)
	{
		ifs->close();
		delete ifs;
	}
}

inline double SRdE(double h, double len, double E)
{
	static const double Cg = 8.85e-05 / twoPi;
	return Cg * pow(E, 4) * h * h * len;
}

AcceleratorModel* MADInterface::ConstructModel()
{
	unique_ptr<DataTable> MADinput(DataTableReaderTFS(filename).Read());

	if(modelconstr != nullptr && appendFlag == false)
	{
		delete modelconstr;
	}
	if(modelconstr == nullptr && appendFlag == false)
	{
		modelconstr = new AcceleratorModelConstructor();
	}

	TypeFactory* factory = new TypeFactory();
	double brho = energy / eV / SpeedOfLight;

	//Loop over all components
	for(size_t i = 0; i < MADinput->Length(); ++i)
	{
		string type = MADinput->Get_s("KEYWORD", i);
		double length = MADinput->Get_d("L", i);

		if(length == 0 && zeroLengths.find(type) != zeroLengths.end())
		{
			MerlinIO::warning() << "Ignoring zero length " << type << ": " << MADinput->Get_s("NAME", i) << endl;
			continue;
		}
		TypeOverrides(MADinput, i);

		if(type == "LINE")
		{
			if(!flatLattice)
			{
				const string& name = MADinput->Get_s("NAME", i);
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
			modelconstr->AppendComponentFrame(ConstructSrot(length, MADinput->Get_s("NAME", i)));
			continue;
		}

		//Determine multipole type by parameters
		vector<AcceleratorComponent*> components = factory->GetInstance(MADinput, energy, brho, i);

		if(inc_sr && (type == "SBEND" || type == "RBEND"))
		{
			energy -= SRdE(MADinput->Get_d("ANGLE", i) / length, length, energy);
			brho = energy / eV / SpeedOfLight;
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
			*log << endl << endl << "final energy = " << energy << " GeV" << endl;
		}
	}

	AcceleratorModel* theModel = modelconstr->GetModel();
	delete modelconstr;
	delete factory;
	modelconstr = nullptr;
	return theModel;
}

void MADInterface::TypeOverrides(unique_ptr<DataTable>& MADinput, size_t index)
{
	string keyword = MADinput->Get_s("KEYWORD", index);
	if(driftTypes.find(keyword) != driftTypes.end())
		MADinput->Set_s("KEYWORD", index, "DRIFT");
	if(keyword == "LCAV")
		MADinput->Set_s("KEYWORD", index, "RFCAVITY");
	if(keyword == "RCOLLIMATOR" || keyword == "ECOLLIMATOR")
		MADinput->Set_s("KEYWORD", index, "COLLIMATOR");
	if(keyword == "RBEND" && MADinput->Get_d("K0L", index))
		MADinput->Set_s("KEYWORD", index, "SBEND");
	if(single_cell_rf && MADinput->Get_s("KEYWORD", index) == "RFCAVITY")
		MADinput->Set_s("KEYWORD", index, "RFCAVITY_SingleCell");
}

string MADInterface::GetMutipoleType(unique_ptr<DataTable>& MADinput, size_t index)
{
	if(!MADinput->Get_d("K0L", index))
		return "SBEND";
	if(!MADinput->Get_d("K1L", index))
		return "QUADRUPOLE";
	if(!MADinput->Get_d("K2L", index))
		return "SEXTUPOLE";
	if(!MADinput->Get_d("K3L", index))
		return "OCTUPOLE";
	if(!MADinput->Get_d("K4L", index))
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
	if(ifs)
	{
		delete ifs;
	}

	ifs = new ifstream(fname.c_str());

	if(!(*ifs))
	{
		MERLIN_ERR << "ERROR opening file " << fname << endl;
		delete ifs;
		abort();
	}

	if(modelconstr == nullptr)
	{
		// first file
		modelconstr = new AcceleratorModelConstructor();
	}

	energy = Pref;

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
			*log << endl << endl << "final energy = " << energy << " GeV" << endl;
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

vector<AcceleratorComponent*> DriftComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	if(length != 0)
		return {new Drift(name, length)};
	else
		return {};
}

vector<AcceleratorComponent*> RBendComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double angle = MADinput->Get_d("ANGLE", id);
	double k1l = MADinput->Get_d("K1L", id);
	double tilt = MADinput->Get_d("TILT", id);
	double h = angle / length;

	SectorBend* bend = new SectorBend(name, length, h, brho * h);

	if(k1l)
		bend->SetB1(brho * k1l / length);

	double e1 = MADinput->Get_d("E1", id);
	double e2 = MADinput->Get_d("E2", id);

	if(e1 != 0 || e2 != 0)
	{
		SectorBend::PoleFace* pf1 = e1 != 0 ? new SectorBend::PoleFace(e1) : nullptr;
		SectorBend::PoleFace* pf2 = e2 != 0 ? new SectorBend::PoleFace(e2) : nullptr;
		bend->SetPoleFaceInfo(pf1, pf2);
	}
	bend->GetGeometry().SetTilt(tilt);

	return {bend};
}

vector<AcceleratorComponent*> SBendComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double angle = MADinput->Get_d("ANGLE", id);
	double k1l = MADinput->Get_d("K1L", id);
	double tilt = MADinput->Get_d("TILT", id);
	double h = angle / length;

	SectorBend* bend = new SectorBend(name, length, h, brho * h);

	if(k1l)
		bend->SetB1(brho * k1l / length);

	double e1 = MADinput->Get_d("E1", id);
	double e2 = MADinput->Get_d("E2", id);

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

vector<AcceleratorComponent*> QuadrupoleComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double k1l = MADinput->Get_d("K1L", id);

	return {new Quadrupole(name, length, brho * k1l / length)};
}

vector<AcceleratorComponent*> SkewQuadrupoleComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy,
	double brho, size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double k1l = MADinput->Get_d("K1L", id);

	return {new SkewQuadrupole(name, length, brho * k1l / length)};
}

vector<AcceleratorComponent*> SextupoleComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double k2l = MADinput->Get_d("K2L", id);

	return {new Sextupole(name, length, brho * k2l / length)};
}

vector<AcceleratorComponent*> SkewSextupoleComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double k2l = MADinput->Get_d("K2L", id);

	return {new SkewSextupole(name, length, brho * k2l / length)};
}

vector<AcceleratorComponent*> OctupoleComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double k3l = MADinput->Get_d("K3L", id);

	return {new Octupole(name, length, brho * k3l / length)};
}

vector<AcceleratorComponent*> YCorComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new YCor(name, length)};
}

vector<AcceleratorComponent*> XCorComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new XCor(name, length)};
}

vector<AcceleratorComponent*> VKickerComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double kick = MADinput->Get_d("VKICK", id);
	double scale;
	if(length > 0)
		scale = brho / length;
	else
		scale = brho;

	return {new YCor(name, length, scale * kick)};
}

vector<AcceleratorComponent*> HKickerComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double kick = MADinput->Get_d("HKICK", id);
	double scale;
	if(length > 0)
		scale = brho / length;
	else
		scale = brho;

	return {new XCor(name, length, -scale * kick)};
}

vector<AcceleratorComponent*> SolenoidComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	double ks = MADinput->Get_d("KS", id);

	return {new Solenoid(name, length, brho * ks / length)};
}

vector<AcceleratorComponent*> RFCavityComponentSingleCell::GetInstance(unique_ptr<DataTable>& MADinput, double energy,
	double brho, size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	// Here we assume an SW cavity
	double freq = MADinput->Get_d("FREQ", id);
	double phase = MADinput->Get_d("LAG", id);
	double volts = MADinput->Get_d("VOLT", id);
	// standing wave cavities need an exact integer of half-wavelengths
	freq *= MHz;
	double lambdaOver2 = SpeedOfLight / freq / 2;
	int ncells = Round(length / lambdaOver2);
	double len1 = ncells * lambdaOver2;

	// adjust phase for cosine-like field
	phase = twoPi * (phase - 0.25);

	double rfcav_len = lambdaOver2;
	double drift_len = length - lambdaOver2;
	ncells = Round(rfcav_len / lambdaOver2);

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

vector<AcceleratorComponent*> RFCavityComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);
	// Here we assume an SW cavity
	double freq = MADinput->Get_d("FREQ", id);
	double phase = MADinput->Get_d("LAG", id);
	double volts = MADinput->Get_d("VOLT", id);
	// standing wave cavities need an exact integer of half-wavelengths
	freq *= MHz;
	double lambdaOver2 = SpeedOfLight / freq / 2;
	int ncells = Round(length / lambdaOver2);
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

vector<AcceleratorComponent*> CrabMarkerComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t id)
{
	double mux = MADinput->Get_d("MUX", id);
	double muy = MADinput->Get_d("MUY", id);
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new CrabMarker(name, length, mux, muy)};
}

vector<AcceleratorComponent*> CrabRFComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new TransverseRFStructure(name, length, 0, 0)};
}

vector<AcceleratorComponent*> CollimatorComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double
	brho, size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new Collimator(name, length)};
}

vector<AcceleratorComponent*> HELComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	return {new HollowElectronLens(name, length, 0, 0, 0, 0, 0)};
}

vector<AcceleratorComponent*> MonitorComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);
	double length = MADinput->Get_d("L", id);

	if(name.substr(0, 2) == "WS")
		return {new RMSProfileMonitor(name, length)};
	else
		return {new BPM(name, length)};
}

vector<AcceleratorComponent*> MarkerComponent::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t
	id)
{
	const string& name = MADinput->Get_s("NAME", id);

	return {new Marker(name)};
}

vector<AcceleratorComponent*> TypeFactory::GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho,
	size_t id)
{
	string type = MADinput->Get_s("KEYWORD", id);
	map<string, getTypeFunc>::iterator itr = componentTypes.find(type);
	if(itr != componentTypes.end())
	{
		return (*itr->second)(MADinput, energy, brho, id);
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
