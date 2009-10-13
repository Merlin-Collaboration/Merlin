/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "merlin_config.h"
#include <algorithm>
// ElementRepository
#include "AcceleratorModel/Implementation/ElementRepository.h"
// deleters
#include "stdext/deleters.h"
// StringPattern
#include "utility/StringPattern.h"

using namespace std;

namespace {

struct MatchID {
public:

    MatchID(const string& idpat,bool negate=false) : pattern(idpat) {}
    bool operator()(const ModelElement* elmnt) const {
        bool rv = pattern(elmnt->GetQualifiedName());
        return neg ? !rv : rv;
    }
private:
    bool neg;
    StringPattern pattern;
};

}; // end anonymous namespace


template class std::set< ModelElement* >;

ElementRepository::~ElementRepository ()
{
    for_each(theElements.begin(),theElements.end(),deleter<ModelElement>());
}

bool ElementRepository::Add (ModelElement* anElement)
{
    pair<iterator,bool> rv = theElements.insert(anElement);
    return rv.second;
}

size_t ElementRepository::Count (const std::string& id) const
{
    return count_if(theElements.begin(),theElements.end(),MatchID(id));
}

size_t ElementRepository::Find (const std::string& id, std::vector<ModelElement*>& elements)
{
    remove_copy_if(theElements.begin(),theElements.end(),
                   back_inserter(elements),MatchID(id,true));
    return elements.size();
}

