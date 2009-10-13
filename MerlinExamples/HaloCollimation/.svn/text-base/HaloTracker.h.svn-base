// class HaloTracking
// ----------------------------------------------------------------------

#ifndef _h_HaloTracker
#define _h_HaloTracker 1

#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/BeamData.h"

class HaloTracker {
public:

    // Construct takes the beamline to track and the initial beam data.
	HaloTracker(AcceleratorModel::Beamline beamline, const BeamData& beamdat);

	// set the halo limits (meter, radian)
	void SetHaloLimits(double x, double xp, double y, double yp, double dp);

	// set the halo limits normalised to the nominal sigma for the transverse
        // coordinates. The momentum cut-off is +-dp as in SetHaloLimits().
	void SetHaloLimitsN(double nx, double nxp, double ny, double nyp, double dp);

	// use STRUCT-type 1/r halo distribution
	void SetSTRUCTlimits(double nx1, double nx2, double ny1, double ny2, double dpp);

    // Run with npart particles in halo
	void Run(size_t npart);

public: // flags

	bool collimate_halo;
	bool dump_particles;
	bool scatter_at_spoiler;

private:

	bool use_struct_dist;
	double n1x,n2x,n1y,n2y;

	double x,xp,y,yp,dp;
	BeamData bdat;
	AcceleratorModel::Beamline bline;
};


#endif
