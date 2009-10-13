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

#include <sstream>
#include "IO/MerlinIO.h"

// MADKeyMap
#include "MADInterface/MADKeyMap.h"
using namespace std;

MADKeyMap::MADKeyMap (const std::string& hstr)
        : has_type(false)
{
    istringstream is(hstr);
    size_t n=0;
    string s;
    while(is>>s) {
        if(s=="TYPE") {
            has_type=true;
        }
        else
            kmap[s]=n++;
    }

#ifndef NDEBUG
    cout<<n<<" column headings identified"<<endl;
#endif

    vals = vector<double>(n,0.0);
}

double MADKeyMap::GetParameter (const std::string& key, bool warn)
{
    key_map::iterator p = kmap.find(key);
    if(p!=kmap.end())
        return vals[p->second];
    else {
        if(warn)
            MerlinIO::warning()<<key<<" not in optics listing. Defaulted to zero"<<endl;
        return 0;
    }
}

void MADKeyMap::ReadRow (std::istream& is)
{
    for(size_t i=0; i<vals.size(); i++)
        is>>vals[i];
}

