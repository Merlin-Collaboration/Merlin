/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TransportMatrix_h
#define TransportMatrix_h 1

#include "merlin_config.h"
#include "TLAS/LinearAlgebra.h"

//	Utility routines for constructing transport matrices (R
//	matrix).

class TransportMatrix
{
public:

    static RealMatrix Drift (double length, RealMatrix& R);
    static RealMatrix& SectorBend (double length, double h, double Kx, RealMatrix& R);
    static RealMatrix& SectorBendT (double length, double h, double Kx, RealMatrix& T);
    static RealMatrix& Quadrupole (double length, double Kx, RealMatrix& R);

    //	Returns in R the 4x4 matrix representing the linear
    //	quadrupole map.
    static RealMatrix& QuadrupoleR (double length, double Kx, RealMatrix& R);

    //	Returns in T the 4x4 matrix representing the
    //	second-order quadrupole map.
    static RealMatrix& QuadrupoleT (double length, double Kx, RealMatrix& T);

    static RealMatrix& Srot (double phi, RealMatrix& R);
    static RealMatrix& Srot (double cosphi, double sinphi, RealMatrix& R);
    static RealMatrix& SrotR4 (double cosphi, double sinphi, RealMatrix& R);
    static RealMatrix& SrotR4 (double phi, RealMatrix& R);
    static RealMatrix& PoleFaceRot (double h, double theta, double fint, double hgap, RealMatrix& R);

    //	Returns a linear matrix for a solenoid having a
    //	longitudinal field strength K0 (1/meter). A solenoid can
    //	also have an optional quadrupole component K1
    //	(1/meter^2). If either entrFringe or exitFringe are
    //	true, then the map contains the effects of the entrance
    //	or exit fringe fields respectively.
    static RealMatrix& Solenoid (double length, double K0, double K1, bool entrFringe, bool exitFringe, RealMatrix& R);

    //	Travelling wave RF cavity (SLAC-91 rev. 2, page 90. )
    static RealMatrix& TWRFCavity (double length, double g, double f, double phi, double E0, bool inc_end_fields, RealMatrix& R);

    //	Standing wave RF cavity
    static RealMatrix& SWRFCavity (int ncells, double g, double f, double phi, double E0, RealMatrix& R);

    //	Returns a linear map and a six-vector (dX) which
    //	approximates an arbitrary 3-D euclidean transformation.
    //	Note the map is in general not symplectic: in the case
    //	of large x- or y-axis rotations, errors may occur.
    static RealMatrix& Arb3DXForm (const RealMatrix& R3d, const RealVector& X, RealMatrix& R, RealVector& dX);
};

inline RealMatrix& TransportMatrix::Quadrupole (double length, double Kx, RealMatrix& R)
{
    return SectorBend(length,0.0,Kx,R);
}

inline RealMatrix& TransportMatrix::Srot (double phi, RealMatrix& R)
{
    return phi==0 ? R : Srot(cos(phi),sin(phi),R);
}

inline RealMatrix& TransportMatrix::SrotR4 (double phi, RealMatrix& R)
{
    return phi==0 ? R : SrotR4(cos(phi),sin(phi),R);
}

#endif
