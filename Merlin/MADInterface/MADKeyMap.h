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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef MADKeyMap_h
#define MADKeyMap_h 1

#include "merlin_config.h"
#include <string>
#include <vector>
#include <iostream>
#include <map>

//      Implementation class for mapping column keys in optics
//      listing to element types during construction.

class MADKeyMap
{
public:

    typedef std::map< std::string , size_t  > key_map;
    struct bad_key {};

    // Constructs a key map from the optics listing line
    // containing the column headings.
    MADKeyMap (const std::string& hstr);

    // Returns the value of the parameter for a specificed key.
    // Throws BAD_KEY if not present.
    virtual double GetParameter (const std::string& key, bool warn=true);

    // Reads in the values for the next row.
    virtual void ReadRow (std::istream& is);

    bool has_type;

 private:

    std::vector<double> vals;
    std::string type_str;
    key_map kmap;
};

#endif
