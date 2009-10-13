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

#ifndef ElementRepository_h
#define ElementRepository_h 1

#include "merlin_config.h"
#include <vector>
#include <string>
#include <set>
// ModelElement
#include "AcceleratorModel/ModelElement.h"

//	Used to store and access all the ModelElement objects
//	associated (contained) by an AcceleratorModel. Primary
//	functions are fast keyed access to ModelElements, and
//	memory management.

class ElementRepository
{
public:
    //	Used to map the element ID's to ModelElement objects in
    //	the repository.

    typedef std::set< ModelElement* > ElementSet;
    typedef ElementSet::iterator iterator;
    typedef ElementSet::const_iterator const_iterator;

    ~ElementRepository ();

    //	Adds the element to the repository. Returns true if
    //	successful, or false if the element already exists in
    //	the repository.
    bool Add (ModelElement* anElement);

    //	Counts the number of elements in the repository with
    //	identifiers which match id. id takes the form
    //	"type.name", where type, name or both can be patterns.
    size_t Count (const std::string& id) const;

    //	Finds and returns in elements all ModelElement objects
    //	whose identifiers match id (see Count() for details of
    //	id).
    size_t Find (const std::string& id, std::vector<ModelElement*>& elements);

    //	Returns the number of ModelElement objects in the
    //	repository.
    size_t Size () const;

    ElementRepository::iterator begin ();
    ElementRepository::iterator end ();
    ElementRepository::const_iterator begin () const;
    ElementRepository::const_iterator end () const;

    ElementSet theElements;
};

inline size_t ElementRepository::Size () const
{
    return theElements.size();
}

inline ElementRepository::iterator ElementRepository::begin ()
{
    return theElements.begin();
}

inline ElementRepository::iterator ElementRepository::end ()
{
    return theElements.end();
}

inline ElementRepository::const_iterator ElementRepository::begin () const
{
    return theElements.begin();
}

inline ElementRepository::const_iterator ElementRepository::end () const
{
    return theElements.end();
}

#endif
