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

#ifndef MerlinIO_h
#define MerlinIO_h 1

#include "merlin_config.h"
#include <iostream>

// MACRO interface
#define MERLIN_OUT  MerlinIO::out()
#define MERLIN_IN   MerlinIO::in()
#define MERLIN_ERR  MerlinIO::error()
#define MERLIN_WARN MerlinIO::warning()

//	Interface for the standard Merlin IO streams.

class MerlinIO
{
public:

    static std::istream& in ();
    static std::ostream& out ();
    static std::ostream& error ();
    static std::ostream& warning ();

    static std::istream* std_in;
    static std::ostream* std_out;
    static std::ostream* std_err;
    static std::ostream* std_warn;
};

inline std::istream& MerlinIO::in ()
{
    return *std_in;
}

inline std::ostream& MerlinIO::out ()
{
    return *std_out;
}

inline std::ostream& MerlinIO::error ()
{
    return *std_err;
}

inline std::ostream& MerlinIO::warning ()
{
    return *std_warn;
}

#endif
