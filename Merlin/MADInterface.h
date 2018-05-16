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
#include <set>

#include "AcceleratorModel.h"

class AcceleratorModelConstructor;
class MADKeyMap;

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
	MADInterface (const std::string& madFileName="", double P0=0);
	~MADInterface();
	/**
	*   Causes the construction of an AcceleratorModel object
	*   based on the MAD optics file.
	*/
	AcceleratorModel* ConstructModel ();

	/**
	*   Sets the log file stream to os.
	*/
	void SetLogFile (ostream& os);

	/**
	*   Turns logging on.
	*/
	void SetLoggingOn ();

	/**
	*   Turns logging off.
	*/
	void SetLoggingOff ();

	/**
	*   If true, all RFCavities will be forced to a length of
	*   wavelength/2 + a Drift of remaining length (LHC MAD tfs table
	*	 bugfix).
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
	void HonourMadStructure (bool flg);

	/**
	*   If true, a flat lattice model in constructed, with no
	*   nested frames.
	*/
	void ConstructFlatLattice (bool flg);

	void ConstructApertures (bool flg);

	/**
	*   Components of type madType are ignored during
	*   construction if their length is zero.
	*/
	void IgnoreZeroLengthType (const string& madType);

	/**
	*   If scaleSR == true, then the magnetic fields of the
	*   magnets are scaled to compensate beam energy losses due
	*   to synchrotron radiation (default = false.) Note that in
	*   this case, the beam energy is the initial energy.
	*/
	void ScaleForSynchRad (bool scaleSR);

	/**
	*   Treats the mad type typestr as a drift.
	*/
	void TreatTypeAsDrift (const std::string& typestr);

	/**
	* Functions for constructing a model from several files.
	* Repeated calls to AppendModel(fname,p) constructs a single
	* model (beamline) from the respective files. The final
	* model is returned using GetModel().
	*/
	void AppendModel(const std::string& fname, double pref);
	AcceleratorModel* GetModel();

	void ConstructNewFrame (const string& name);
	void EndFrame (const string& name);

protected:
	double energy;
	std::string filename;
	ifstream *ifs;
	ostream* log;

	bool logFlag;
	bool flatLattice;
	bool honMadStructs;
	bool incApertures;
	bool inc_sr;

	std::set<std::string> zeroLengths;
	std::set<std::string> driftTypes;

	AcceleratorModelConstructor* ctor;
	MADKeyMap* prmMap;

	double ReadComponent ();
	void Initialise();
	double z;	///Distance along the lattice

	bool single_cell_rf;
};

inline void MADInterface::SetLogFile (ostream& os)
{
	log=&os;
}

inline void MADInterface::SetLoggingOn ()
{
	logFlag=true;
}

inline void MADInterface::SetLoggingOff ()
{
	logFlag=false;
}

inline void MADInterface::HonourMadStructure (bool flg)
{
	honMadStructs=flg;
}

inline void MADInterface::ConstructFlatLattice (bool flg)
{
	flatLattice=flg;
}

inline void MADInterface::ConstructApertures (bool flg)
{
	incApertures = flg;
}

inline void MADInterface::ScaleForSynchRad (bool scaleSR)
{
	inc_sr = scaleSR;
}
#endif
