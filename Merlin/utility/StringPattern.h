/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/01 12:34:18 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef StringPattern_h
#define StringPattern_h 1

#include "merlin_config.h"
#include <string>
#include <vector>
#include <utility>
#include <iostream>

using std::string;
using std::ostream;

//	A simple string pattern which represents a reduced set
//	of the standard UNIX regular expression syntax.
//	Currently supported is the single wild-card '*', which
//	matches any number (including zero) of any character.
//  StringPatter also supports OR: pat1|pat2|pat3 etc.

class StringPattern
{
public:

    //	Constructor taking the string pattern.
    StringPattern (const std::string& s);

    //	Match operation. Returns true if the pattern matched the
    //	string.
    bool Match (const std::string& s) const;

    //	Operator form of match().
    bool operator () (const std::string& s) const;

    //	Outputs to os the original pattern.
    friend ostream& operator << (ostream& os, const StringPattern& pattern);

    operator string () const;

    bool operator < (const StringPattern& rhs) const;
    bool operator == (const StringPattern& rhs) const;

    //	The wildcard character (default ='*').
    static char wcchar;

private:

    //	A vector of literal strings. The pattern is split into
    //	literals by the presence of a '*' character.
    std::vector<string> patterns;

	// list of sub-patters for an OR'd pattern:
	std::vector<StringPattern> orpatterns;

    //	Used to indicate if the first/last character of the
    //	string pattern were wildcards.
    std::pair<bool,bool> wcterm;

    //	Set to true if this pattern is a simple literal string
    //	(no wild cards).
    bool isLiteral;

    std::string str;
};

inline bool StringPattern::operator () (const std::string& s) const
{
    return Match(s);
}

inline StringPattern::operator string () const
{
    return str;
}

inline bool StringPattern::operator < (const StringPattern& rhs) const
{
    return str<rhs.str;
}

inline bool StringPattern::operator == (const StringPattern& rhs) const
{
    return str==rhs.str;
}

// Utility function for splitting a single string into a list of delimited strings.
// The return pair indicates if the first and last characters of str are dchar.
std::pair<bool,bool> DelimitString(const string& str, char dchar, std::vector<string>& result);

#endif
