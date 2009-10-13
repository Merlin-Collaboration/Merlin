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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef SimpleATL_h
#define SimpleATL_h 1

#include "merlin_config.h"
#include <iostream>
#include <vector>
// AcceleratorSupport
#include "AcceleratorModel/Supports/AcceleratorSupport.h"

class RandGenerator;

//	Represents a simple ATL model of ground motion. On each
//	time step dT, SimpleATL applies a random vertical (and
//	horizontal if included) displacement to a set of
//	AcceleratorSupports. The variance of the motion is
//
//	                 v = A.dT.L,
//
//	where A is the ATL constant, dT the time step and L the
//	distance between supports.
//
//	Note that at present only the arc distance between
//	successive supports  is used. For curved geometries this
//	may introduce an error.

class SimpleATL
{
public:
    //	Constructor taking the A constant and the list of
    //	support structures.
    SimpleATL (double anA, const AcceleratorSupportList& supports, double vrms=0);

    ~SimpleATL ();

    //	Reset the ground motion to zero Note this resets the
    //	offset of all the AcceleratorSupports, and resets the
    //	internal clock to zero.
    void Reset ();

    //	Perform a single step of dt seconds. Returns the current
    //	simulated time.
    double DoStep (double dt);

    //	Record the (x,y,z) offset of all the supports to the
    //	specified stream.
    void RecordOffsets (std::ostream& os) const;

    //	Returns the current simulated time (in seconds).
    double GetTime () const;

    //	Sets the random seed to nseed.
    void SetRandomSeed (unsigned int nseed);

    //	Returns the current random seed
    unsigned int GetRandomSeed () const;

    //	Resets the random generator with the current random seed.
    void ResetRandomSeed ();

private:

    double t;
    double A;
    unsigned seed;

    // Uncorrelated white-noise vibration variance
    double vv;

    // vector to store ATL ground motion
    std::vector<double> atlgm;

    AcceleratorSupportList theSupports;

    RandGenerator* rg;
};

#endif
