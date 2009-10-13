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

#ifndef ATL2D_h
#define ATL2D_h 1

#include "merlin_config.h"
#include <iostream>
#include <vector>
// AcceleratorSupport
#include "AcceleratorModel/Supports/AcceleratorSupport.h"
// TLASimp
#include "TLAS/LinearAlgebra.h"

class RandGenerator;

//	Represents a 2D ATL model of ground motion. On each
//	time step dT, ATL2D applies a random vertical (and
//	horizontal if included) displacement to a set of
//	AcceleratorSupports. The variance of the motion is
//
//	                 v = A.dT.L,
//
//	where A is the ATL constant, dT the time step and L the
//	distance between supports.
//
//  ATL2D calculates correctly the correlation between the support point
//  on a 2D plane, such that A.dT.L  holds for any two points, where L is
//  the direct distance between those two points.
//
//

class ATL2D
{
public:

    enum ATLMode {increment, absolute};

    //	Constructor taking the A constant and the list of
    //	support structures.
    ATL2D (double anA, const AcceleratorSupportList& supports, const Point2D refPoint=Point2D(0,0), ifstream* evecTFile=0, ifstream* evalFile=0);

    ~ATL2D ();

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

    bool SetATLMode (const ATLMode mode);

    bool SetVibration (const double vrms);

    void RecordEigenSystem (ofstream* evecTFile, ofstream* evalFile);

private:

    double t;
    double A;
    unsigned seed;

    // Uncorrelated white-noise vibration variance
    double vv;

    AcceleratorSupportList theSupports;
    RandGenerator* rg;

    RealMatrix evecsT;
    RealVector evals;
    ATLMode atlMode;

    double Distance(const int n1, const int n2);
    double Distance(const int n1, const Point2D x2);
};

#endif
