/////////////////////////////////////////////////////////////////////////
// Class ResponseMatrixGenerator
// Utility class for constructing the response matrix for a given set of
// correctors and BPMs. Once constructed, the response matrix and reference
// trajectory can be calculated for any valid beam state.
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

#ifndef _h_ResponseMatrixGenerator
#define _h_ResponseMatrixGenerator

#include "TLAS/TLAS.h"
#include "Channels/Channels.h"
#include "Accelerator.h"

class ResponseMatrixGenerator {
public:

	ResponseMatrixGenerator(Accelerator* acc, const ROChannelArray& bpms, RWChannelArray& cors, double eps = 1.0e-06);
	
	// Generate the response matrix for the ns-th beam state
	const RealMatrix& Generate(size_t ns);

	const RealMatrix& GetMatrix() const;
	const RealVector& GetReference() const;


private:

	Accelerator* acc;
	RWChannelArray& cors;
	const ROChannelArray& bpms;
	double eps;

	RealVector data0;
	RealMatrix M;
};

#endif