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

#ifndef _h_MatrixPrinter
#define _h_MatrixPrinter

#include "merlin_config.h"
#include <iostream>
#include "TLAS/LinearAlgebra.h"
#include "utility/OPFormat.h"


void MatrixForm(const RealMatrix&, std::ostream&, const OPFormat& =OPFormat());


#endif
