/////////////////////////////////////////////////////////////////////////
// Class Accelerator
// Represents the physical accelerator. Provides the primary interface
// to the tuning application to the underlying accelerator model.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_Accelerator
#define _h_Accelerator

#include <vector>
#include <string>
#include <utility>

#include "CommonDataStructures.h"
#include "BeamModel/Bunch.h"
#include "BeamModel/BeamData.h"

class RWChannelArray;
class ROChannelArray;
class BeamDynamicsModel;
class AcceleratorModel;

// Represents the physical accelerator. Provides the primary interface
// to the tuning application to the underlying accelerator model.

typedef std::vector<size_t> IntegerArray;

class Accelerator {	
public:

	// Indicates the working plane.
	enum Plane {
		x_only,
		y_only,
		x_and_y
	};

	Accelerator(const std::string& name, AcceleratorModel*, BeamData*);
	~Accelerator();

	// Return the name of this accelerator
	const std::string& GetName() const;

	// Set the beam dynamics model to use for tracking
	void SetBeamDynamicsModel(BeamDynamicsModel*);

	// Control incremental tracking
	// Note that setting the flag to true does not guarantee
	// that incremental tracking will be used.
	void AllowIncrementalTracking(bool);

	// Sets the current beamline segment being tuned. This function affects
	// which range of BPMs and Correctors are used by GetMonitorChannels and
	// GetCorrectorChannels. 
	void SetActiveBeamlineSegment(const DFS_Segment& seg);	

	// Track the beam corresponding to state n. 
	// The state index n is used when incremental 
	// tracking is implemented.
	void TrackBeam(size_t n);	

	// Construct a new bunch and track it through the entire
	// model. Cached bunch state and active segment are ignored
	void TrackNewBunchThroughModel();

	// Return in bpmChannels the BPM channels for the specified plane
	// and the active beamline segment. Number of channels is returned.
	size_t GetMonitorChannels(Plane xy, ROChannelArray& bpmChannels);

	// Return in corrChannels the corrector channels for the specified plane
	// and the active beamline segment. Number of channels is returned.
	size_t GetCorrectorChannels(Plane xy, RWChannelArray& corrChannels);

	// Returns in klys the Klystrons for the Accelerator. Number of 
	// klystrons is returned.
	size_t GetKlystrons(KlystronArray& klys);

	// Initialises the tracking engine with the specified number of energy states.
	// Returns in refplist the reference particles for each initial beam for each state.
	void InitialiseTracking(size_t nstates, ReferenceParticleArray& refplist);
	
	// Returns in indecies the beamline indecies of the
	// components matching cpattern. Returns the length of indecies
	// on exit.
	size_t GetBeamlineIndecies(const std::string& cpat, IntegerArray& indecies) const;

	// Returns the beamline indecies for the entire model
	DFS_Segment GetBeamlineRange() const;

	// Set the BPM single-shot resolution
	void SetBPMresolution(double rms);

protected:
	AcceleratorModel* itsAccModel;	

private:
	const std::string itsName;
	BeamDynamicsModel* itsTracker;	
	BeamData* beam0;	

	// Ordered list of bunch objects used for incremental tracking.
	class CachedBunch {
	public:
		size_t location;
		mutable Bunch* bunch;
		explicit CachedBunch(Bunch* aBunch)
			: location(0),bunch(aBunch)
		{}
		CachedBunch(const CachedBunch& rhs)
			: location(rhs.location),bunch(rhs.bunch){
			rhs.bunch=0;
		}

		~CachedBunch() {
			if(bunch)
				delete bunch; 
		}
	};

	std::vector<CachedBunch> cachedBunches;
	DFS_Segment currentSegment;
	bool allowIncrTracking;

	// Copy construction/assignment not allowed.
	Accelerator(const Accelerator&);
	void operator=(const Accelerator&);
};

#endif

