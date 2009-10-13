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

#ifndef _h_RotationMatrix
#define _h_RotationMatrix 1

#include "TLAS/LinearAlgebra.h"

// functions for constructing 3x3 rotation arrays

inline RealMatrix DECL_XROT(double c, double s)
{
    RealMatrix m(IdentityMatrix(3));
    m(1,1)=m(2,2)=c;
    m(1,2)=s;
    m(2,1)=-s;
    return m;
}

inline RealMatrix DECL_YROT(double c, double s)
{
    RealMatrix m(IdentityMatrix(3));
    m(0,0)=m(2,2)=c;
    m(0,2)=-s;
    m(2,0)=s;
    return m;
}

inline RealMatrix DECL_ZROT(double c, double s)
{
    RealMatrix m(IdentityMatrix(3));
    m(0,0)=m(1,1)=c;
    m(0,1)=s;
    m(1,0)=-s;
    return m;
}

#endif // _h_RotationMatrix
