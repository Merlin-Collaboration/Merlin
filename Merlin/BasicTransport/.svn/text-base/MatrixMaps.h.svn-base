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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef MatrixMaps_h
#define MatrixMaps_h 1

#include "merlin_config.h"
#include <cassert>
// PSTypes
#include "BeamModel/PSTypes.h"
// LinearAlgebra
#include "TLAS/LinearAlgebra.h"

//	Base class for linear matrix maps. LinMtrxMap is
//	intended only as a common base for class RMtrx and class
//	RdpMtrx. It cannot be instantiated.

class LinMtrxBase
{
public:

    //	Return the number of degrees of freedom for this matrix
    //	(dimension/2).
    int GetNDF () const;

    //	Return the reference momentum in GeV/c (returns 0 if the
    //	map is energy independent).
    double GetRefMomentum () const;

    //	Returns true if this map is energy independent.
    bool EnergyIndependent () const;

    RealMatrix R;

protected:

    LinMtrxBase (double P, int ndf);

    //	Construction directly from an arbitrary matrix.
    LinMtrxBase (const RealMatrix& RR, double P0);

    //	Returns the momentum error re-scaled to the reference
    //	momentum.
    double scaledp (double P, double dp) const;

    double Pref;
};

//	RMtrx represents a linear map which is implemented by a
//	matrix (R matrix). The map can act on single vectors or
//	first- and second-order moments, or arrays of the above.
//	The dimensionality of the map is specified in the number
//	of degrees of freedom (1,2 or 3).  In addition, RMtrx
//	objects are associated with a specific reference
//	momentum, about which the map (R) is calculated.
//	Functions are provided for mapping vectors or moments
//	whose reference momenta are different from the RMtrx
//	reference. RMtrx objects can be multiplied together and
//	inverted.

class RMtrx : public LinMtrxBase
{
public:

    //	Constructor taking the reference momentum and the number
    //	of degrees of freedom (ndf) for the matrix. The  matrix
    //	is dimensioned but not initialised.
    explicit RMtrx (int ndf = 3, double P = 0);

    //	Construction from an aribrary matrix. The matrix should
    //	be square with dimension 2n, where n is 1,2 or 3.
    RMtrx (const RealMatrix& RR, double P0 = 0);

    // Set the reference momentum
    void SetRefMomentum(double p) { Pref = p; }

    //	Transforms the vector x by R such that x->R.x. The
    //	second form explicitly specifies the reference momentum
    //	(GeV/c) for the vector. Returns x.
    PSvector& Apply (PSvector& x) const;
    PSvector& Apply (PSvector& x, double p0) const;

    //	Transforms each vector in xa by R such that x->R.x. The
    //	second form explicitly specifies the reference momentum
    //	(GeV/c) for the vectors. Returns xa.
    PSvectorArray& Apply (PSvectorArray& xa) const;
    PSvectorArray& Apply (PSvectorArray& xa, double p0) const;

    //	Transform  sigma by R. If X,S represent the first- and
    //	second-order moments respectively, then X->R.X and
    //	S->R.S.R'. The second form explicitly specifies the
    //	reference momentum (GeV/c)  for the moments. Returns
    //	sigma.
    PSmoments& Apply (PSmoments& sigma) const;
    PSmoments& Apply (PSmoments& sigma, double p0) const;

    //	Transforms the moments in sigmaArray by R (see
    //	apply(PSmoments&) for details). The second form
    //	explicitly specifies the reference momentum (GeV/c) for
    //	the moments. Returns sigmaArray.
    PSmomentsArray& Apply (PSmomentsArray& sigmaArray) const;
    PSmomentsArray& Apply (PSmomentsArray& sigmaArray, double p0) const;

    //	Invert the matrix.
    RMtrx& Invert ();

    //	Dot operation. (*this)*=R represents (*this)->R.(*this).
    RMtrx& operator *= (const RMtrx& rhs);

    // array-like accessor
    double operator()(int i, int j) const { return R(i-1,j-1); }
};

//	An approximate linear map which contains a first-order
//	chromatic correction. The matrix used for tracking is
//	calculated as R0+R1*dp, where R0 is the matrix
//	calculated about the reference momentum, R1 is the first
//	order derivative wrt momentum and dp is the momentum
//	error. For more information see RMtrx.

class RdpMtrx : public LinMtrxBase
{
public:

    //	Constructor taking the reference momentum and the number
    //	of degrees of freedom (ndf) for the matrix. The
    //	matrices are dimensioned but not initialised.
    explicit RdpMtrx (int ndf = 3, double P0 = 0);

    //	Constructor from two arbitrary matrices. Both matrices
    //	should be square and of equal dimension, which should be
    //	2n where n=1,2, or 3.
    RdpMtrx (const RealMatrix& RR, const RealMatrix& TT, double P0 = 0);

    //	Transforms the vector x by R such that x->R.x. The
    //	second form explicitly specifies the reference momentum
    //	(GeV/c) for the vector. Returns x.
    PSvector& Apply (PSvector& x) const;
    PSvector& Apply (PSvector& x, double p0) const;

    //	Transforms each vector in xa by R such that x->R.x. The
    //	second form explicitly specifies the reference momentum
    //	(GeV/c) for the vectors. Returns xa.
    PSvectorArray& Apply (PSvectorArray& xa) const;
    PSvectorArray& Apply (PSvectorArray& xa, double p0) const;

    //	Transform  sigma by R. If X,S represent the first- and
    //	second-order moments respectively, then X->R.X and
    //	S->R.S.R'. The second form explicitly specifies the
    //	reference momentum (GeV/c)  for the moments. Returns
    //	sigma.
    PSmoments& Apply (PSmoments& sigma) const;
    PSmoments& Apply (PSmoments& sigma, double p0) const;

    //	Transforms the moments in sigmaArray by R (see
    //	apply(PSmoments&) for details). The second form
    //	explicitly specifies the reference momentum (GeV/c) for
    //	the moments. Returns sigmaArray.
    PSmomentsArray& Apply (PSmomentsArray& sigmaArray) const;
    PSmomentsArray& Apply (PSmomentsArray& sigmaArray, double p0) const;

    //	Second-order momemtum derivative matrix (T matrix).
    RealMatrix T;
};

inline LinMtrxBase::LinMtrxBase (double P, int ndf)
        : R(2*ndf,2*ndf),Pref(P)
{}

inline int LinMtrxBase::GetNDF () const
{
    return R.nrows()/2;
}

inline double LinMtrxBase::GetRefMomentum () const
{
    return Pref;
}

inline bool LinMtrxBase::EnergyIndependent () const
{
    return Pref==0;
}

inline double LinMtrxBase::scaledp (double P, double dp) const
{
    return P*(1+dp)/Pref-1.0;
}

inline RMtrx::RMtrx (int ndf, double P)
        : LinMtrxBase(P,ndf)
{}

inline RMtrx::RMtrx (const RealMatrix& RR, double P0)
        : LinMtrxBase(RR,P0)
{}

inline RdpMtrx::RdpMtrx (int ndf, double P0)
        : LinMtrxBase(P0,ndf),T(2*ndf,2*ndf)
{}

inline RdpMtrx::RdpMtrx (const RealMatrix& RR, const RealMatrix& TT, double P0)
        : LinMtrxBase(RR,P0),T(TT)
{
    assert(TT.nrows()==RR.nrows() && TT.ncols()==RR.ncols());
}

#endif
