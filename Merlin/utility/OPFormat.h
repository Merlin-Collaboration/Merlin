/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef OPFormat_h
#define OPFormat_h 1

#include "merlin_config.h"
#include <iostream>
#include <string>

using std::string;
using std::ios_base;
using std::ostream;

//	Floating Point Format. Utility class for formatting and
//	outputting floating point numbers to an ostream, without
//	modifying the state of the ostream.

class OPFormat
{
public:
    
    explicit OPFormat (int p = 6)
        :prc(p),wdt(0),fmt(),adjust(ios_base::right),allowovf(true)
    {}
    
    OPFormat& scientific ()
    {
        fmt=ios_base::scientific;
        return *this;
    }
    
    OPFormat& fixed ()
    {
        fmt=ios_base::fixed;
        return *this;
    }
    
    OPFormat& general ()
    {
        fmt &= ~(ios_base::scientific & ios_base::fixed);
        return *this;
    }
    
    OPFormat& precision (int p)
    {
        prc=p;
        return *this;
    }
    
    OPFormat& width (int w)
    {
        wdt=w;
        return *this;
    }
    
    OPFormat& left ()
    {
        adjust=ios_base::left;
        return *this;
    }
    
    OPFormat& right ()
    {
        adjust=ios_base::right;
        return *this;
    }
    
    OPFormat& internal ()
    {
        adjust=ios_base::internal;
        return *this;
    }
    
    OPFormat& overflow (bool flg)
    {
        allowovf=flg;
        return *this;
    }
    
    string operator () (double val) const;
    string operator () (int val) const;
    string operator () (const string& val) const;
    
    friend ostream& operator << (ostream& os, const OPFormat& fmt);
    
    int prc;
    int wdt;
    ios_base::fmtflags fmt;
    ios_base::fmtflags adjust;
    bool allowovf;
};

#endif
