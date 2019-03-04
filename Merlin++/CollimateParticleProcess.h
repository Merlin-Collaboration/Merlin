/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CollimateParticleProcess_h
#define CollimateParticleProcess_h 1

#include <map>
#include <set>
#include <list>
#include <vector>

#include "merlin_config.h"
#include "ParticleBunchProcess.h"
#include "PSTypes.h"
#include "CollimationOutput.h"
#include "MerlinException.h"

#define COLL_AT_ENTRANCE 1
#define COLL_AT_CENTER 2
#define COLL_AT_EXIT 4

namespace ParticleTracking
{

class ExcessiveParticleLoss: public MerlinException
{
public:
	ExcessiveParticleLoss(const string& c_id, double threshold, size_t nlost, size_t nstart);
};

/**
 * A process which effectively collimates the particle
 * beam, given the aperture of the current  component.
 * Collimation within a component can occur at one or more
 * positions:

 * COLL_AT_ENTRANCE,
 * COLL_AT_CENTER,
 * COLL_AT_EXIT
 */
class CollimateParticleProcess: public ParticleBunchProcess
{
public:

	/**
	 * Used to generate unique filenames for particle loss data files.
	 */
	typedef map<string, int> IDTBL;

	/**
	 * Constructor taking the collimation mode, and the output
	 * stream pointer to which to print the results. mode can
	 * be a logical OR combination of the collimation modes. A
	 * nullptr for osp (default) suppresses output.
	 */
	CollimateParticleProcess(int priority, int mode, std::ostream* osp = nullptr);

	~CollimateParticleProcess();

	/**
	 * Initialise this process with the specified Bunch. If
	 * bunch is not a ParticleBunch object, the process becomes
	 * inactive.
	 */
	virtual void InitialiseProcess(Bunch& bunch);

	/**
	 * Sets the current accelerator component.
	 */
	virtual void SetCurrentComponent(AcceleratorComponent& component);

	/**
	 * Preform the process for the specified step ds.
	 */
	virtual void DoProcess(double ds);

	/**
	 * Returns the current maximum step length for this process.
	 * @return Current maximum step length
	 */
	virtual double GetMaxAllowedStepSize() const;

	/**
	 * If set to true, the process scatters the particles in
	 * energy and angle at a Collimator element, if the particle is
	 * outside the aperture.
	 */
	void ScatterAtCollimator(bool tf);

	/**
	 * If flg is true, then files are generated containing the
	 * lost (collimated) particles. The file names have the
	 * form {fprefix.}type.id.n.loss, where the string fprefix
	 * is optional, type.id is the qualified name of the
	 * element where the particle loss occurs, and n is an
	 * occurrence count for like-named elements (starting with
	 * n=1).
	 */
	void CreateParticleLossFiles(bool flg, string fprefix = "");

	/**
	 * If index is true, then the initial  particles are
	 * sequential indexed (1..n). These index values for each
	 * particle are then maintained, and output during any
	 * particle output operation (as the first column). The
	 * indexing allows particles to be traced back to the
	 * original coordinates.
	 */
	void IndexParticles(bool index);
	void IndexParticles(std::list<size_t>& anIndex);
	const std::list<size_t>& GetIndexes() const;

	/**
	 * Sets the threshold for particle loss before the process
	 * throws ParticleLossThreshold exception. The value is in
	 * % of the initial particle number (default = 100%).
	 */
	void SetLossThreshold(double losspc);

	/**
	 * Set the log stream for the process. A nullptr
	 * turns logging off.
	 */
	void SetLogStream(std::ostream* anOs);

	/**
	 * Enable collimator jaw imperfections
	 */
	void EnableImperfections(bool);

	virtual double GetOutputBinSize() const;
	virtual void SetOutputBinSize(double);

	virtual void SetCollimationOutput(CollimationOutput* odb)
	{
		CollimationOutputVector.push_back(odb);
		CollimationOutputSet = true;
	}

	std::vector<CollimationOutput*> CollimationOutputVector;
	std::vector<CollimationOutput*>::iterator CollimationOutputIterator;

protected:

	int cmode;
	std::ostream* os;
	bool createLossFiles;
	string file_prefix;
	double lossThreshold;
	size_t nstart;
	std::list<size_t>* pindex;

	IDTBL idtbl;

	/**
	 * The input array
	 */
	PSvectorArray InputArray;

	double s_total;
	double s;
	double next_s;

	/**
	 * physical length
	 */
	double len;

	bool is_collimator;
	bool at_entr;
	bool at_cent;
	bool at_exit;

	size_t nlost;

	/**
	 * 0 when no CollimationOutput is set
	 * @retval 0 No CollimationOutput set
	 */
	bool CollimationOutputSet;

	double GetBinSize() const
	{
		return bin_size;
	}

	int ColParProTurn;
	std::string FirstElementName;
	double FirstElementS;
	bool FirstElementSet;

private:

	virtual void DoCollimation();
	void SetNextS();
	virtual void DoOutput(const PSvectorArray& lostb, const std::list<size_t>& lost_i);
	void bin_lost_output(const PSvectorArray& lostb);

	bool scatter;
	double bin_size;        /// size of bins
	double step_size;       /// step size to handle when the collimator cannot be split equally into bins
	bool Imperfections;

	double Xr; /// radiation length
	virtual bool DoScatter(Particle&);

	/**
	 * A list of particles we want to use in the input array
	 */
	std::vector<unsigned int> LostParticlePositions;
};

inline void CollimateParticleProcess::CreateParticleLossFiles(bool flg, string fprefix)
{
	createLossFiles = flg;
	file_prefix = fprefix;
}

inline void CollimateParticleProcess::SetLogStream(std::ostream* anOs)
{
	os = anOs;
}

inline const std::list<size_t>& CollimateParticleProcess::GetIndexes() const
{
	return *pindex;
}

inline void CollimateParticleProcess::ScatterAtCollimator(bool tf)
{
	scatter = tf;
}

} // end namespace ParticleTracking

#endif
