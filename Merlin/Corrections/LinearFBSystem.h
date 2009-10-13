/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef LinearFBSystem_h
#define LinearFBSystem_h 1

#include "merlin_config.h"

#include "TLAS/TLAS.h"
#include <queue>
#include <list>
// Channels
#include "Channels/Channels.h"
// LinearAlgebra
#include "TLAS/LinearAlgebra.h"

using namespace TLAS;

//	A simple linear feedback correction algorithm. On each
//	application, the actuator channels (A) are incremented
//	using the following linear equation:
//
//	A = A-g*(Mi*(S-S0))
//
//	where g is the gain, S are the current signal values and
//	S0 are the desired signal values (set points). Mi is a
//	speudo-inverse matrix of the response matrix M defined by
//
//	S=M*A
//
//	Mi is calculated using SVD.

class LinearFBSystem
{
public:

    LinearFBSystem (std::vector<ROChannel*>& sigs, std::vector<RWChannel*>& acts, const RealMatrix& M);
    LinearFBSystem (ROChannelArray& sigs, RWChannelArray& acts, const RealMatrix& M);

    ~LinearFBSystem ();

    void SignalsToSetpoints ();
    void StoreActuators () const;
    void RestoreActuators ();
    void SetResponseMatrix (const RealMatrix& M);
    void SetGain (double g);
    double GetGain () const;
    void Apply ();
    void SetSetpoints (const RealVector& S0);
    double GetActuatorRMS () const;
    double GetSignalRMS () const;
    int GetNumSignals () const;
    int GetNumActuators () const;
    void SetPulseDelay(int n);

private:

    double gain;
    ROChannelArray signals;
    RWChannelArray actuators;
    RealVector setpoints;

    LinearFBSystem(const LinearFBSystem &right);
    const LinearFBSystem & operator=(const LinearFBSystem &right);

    mutable RealVector* cached_actuators;

    // to allow for possible actuator pulse delays, we use a queue
    mutable std::queue<RealVector>* actuatorQueue;

    SVDMatrix< double >* Mi;
};

#endif
