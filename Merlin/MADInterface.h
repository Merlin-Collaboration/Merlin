/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MADInterface_h
#define MADInterface_h 1

#include "merlin_config.h"
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <memory>
#include "AcceleratorModel.h"
#include "DataTable.h"

class AcceleratorModelConstructor;

using std::ifstream;
using std::ostream;

/**
 *      Class used to construct a MERLIN model from a MAD optics
 *      output listing. The class now automatically  identifies
 *      the column parameters, and associates them with the
 *      constructed element types. If an element type is defined
 *      for which a required parameter is not present in the
 *      column headings, the parameter is set to zero and a
 *      warning is issued.
 */

class MADInterface
{
public:
	/**
	 *  Constructor taking the name of the MAD optics file, and
	 *  the momentum in GeV/c.
	 */
	MADInterface(const std::string& madFileName = "", double P0 = 0);

	/**
	 *  Destructor
	 */
	~MADInterface();

	AcceleratorComponent* ComponentTypeMap(string& type);

	/**
	 *   Causes the construction of an AcceleratorModel object
	 *   based on the MAD optics file.
	 */
	AcceleratorModel* ConstructModel();

	/**
	 *   Sets the log file stream to os.
	 */
	void SetLogFile(ostream& os);

	/**
	 *   Turns logging on.
	 */
	void SetLoggingOn();

	/**
	 *   Turns logging off.
	 */
	void SetLoggingOff();

	/**
	 * Function to return corresponding multipole string *
	 */
	string GetMutipoleType(unique_ptr<DataTable>& MADinput, size_t id);

	/**
	 * Function to define all type overrides
	 */
	void TypeOverrides(unique_ptr<DataTable>& MADinput, size_t index);

	/**
	 *   If true, all RFCavities will be forced to a length of
	 *   wavelength/2 + a Drift of remaining length (
	 *   //LHC MAD tfs table bugfix!!!
	 */
	void SetSingleCellRF(bool scrf)
	{
		single_cell_rf = scrf;
	}

	/**
	 *   If true, all LINE constructs in the MAD optics output
	 *   are constructed in the model. If false, only those
	 *   prefixed X_, where X is M, S, or G are constructed.
	 */
	void HonourMadStructure(bool flg);

	/**
	 *   If true, a flat lattice model in constructed, with no
	 *   nested frames.
	 */
	void ConstructFlatLattice(bool flg);

	/**
	 *   Components of type madType are ignored during
	 *   construction if their length is zero.
	 */
	void IgnoreZeroLengthType(const string& madType);

	/**
	 *   If scaleSR == true, then the magnetic fields of the
	 *   magnets are scaled to compensate beam energy losses due
	 *   to synchrotron radiation (default = false.) Note that in
	 *   this case, the beam energy is the initial energy.
	 */
	void ScaleForSynchRad(bool scaleSR);

	/**
	 *   Treats the mad type typestr as a drift.
	 */
	void TreatTypeAsDrift(const std::string& typestr);

	void AppendModel(const std::string& fname, double pref);
	AcceleratorModel* GetModel();
	AcceleratorModelConstructor* GetModelConstructor();

	void ConstructNewFrame(const string& name);
	void EndFrame(const string& name);

	double GetEnergy()
	{
		return energy;
	}

	void SetEnergy(double newenergy)
	{
		energy = newenergy;
	}

	bool GetSynchRadFlag()
	{
		return inc_sr;
	}

	double energy;
	bool inc_sr;
	bool flatLattice;
	double z;   ///Distance along the lattice
	bool single_cell_rf;

protected:

	std::string filename;
	ifstream *ifs;
	ostream* log;

	bool logFlag;
	bool honMadStructs;
	bool appendFlag;

	std::set<std::string> zeroLengths;
	std::set<std::string> driftTypes;

	AcceleratorModelConstructor* modelconstr;
	AcceleratorComponent* currentcomponent;
};

inline void MADInterface::SetLogFile(ostream& os)
{
	log = &os;
}

inline void MADInterface::SetLoggingOn()
{
	logFlag = true;
}

inline void MADInterface::SetLoggingOff()
{
	logFlag = false;
}

inline void MADInterface::HonourMadStructure(bool flg)
{
	honMadStructs = flg;
}

inline void MADInterface::ConstructFlatLattice(bool flg)
{
	flatLattice = flg;
}

inline void MADInterface::ScaleForSynchRad(bool scaleSR)
{
	inc_sr = scaleSR;
}

typedef AcceleratorComponent* (*getTypeFunc)(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);

class TypeFactory
{
public:
	static map<string, getTypeFunc> componentTypes;
	AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);

};

class TypeFactoryInit
{
	static TypeFactoryInit init;
public:
	TypeFactoryInit();
};

class DriftComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class RBendComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SBendComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class QuadrupoleComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SkewQuadrupoleComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SextupoleComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SkewSextupoleComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class OctupoleComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class YCorComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class XCorComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class VKickerComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class HKickerComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SolenoidComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class RFCavityComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class CollimatorComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class CrabMarkerComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class CrabRFComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class HELComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class MonitorComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class MarkerComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class LineComponent: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

class SROTComponenet: public AcceleratorComponent
{
public:
	static AcceleratorComponent* GetInstance(unique_ptr<DataTable>& MADinput, double energy, double brho, size_t id);
};

#endif
