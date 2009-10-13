/////////////////////////////////////////////////////////////////////////
// ILCDFS IO support implementation
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
#include "ILCDFS_IO.h"

std::ostream* dfs_trace::os = &std::cout;
dfs_trace::trace_level dfs_trace::verbosity = dfs_trace::level_3;
