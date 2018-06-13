/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "merlin_config.h"
#include <algorithm>
#include "ElementRepository.h"
#include "deleters.h"
#include "StringPattern.h"

using namespace std;

namespace
{

struct MatchID
{
public:

	MatchID(const string& idpat, bool negate = false) :
		pattern(idpat)
	{
	}
	bool operator()(const ModelElement* elmnt) const
	{
		bool rv = pattern(elmnt->GetQualifiedName());
		return neg ? !rv : rv;
	}
private:
	bool neg;
	StringPattern pattern;

};

} // end anonymous namespace

template class std::set<ModelElement*>;

ElementRepository::~ElementRepository()
{
	for_each(theElements.begin(), theElements.end(), deleter<ModelElement>());
}

bool ElementRepository::Add(ModelElement* anElement)
{
	pair<iterator, bool> rv = theElements.insert(anElement);
	return rv.second;
}

size_t ElementRepository::Count(const std::string& id) const
{
	return count_if(theElements.begin(), theElements.end(), MatchID(id));
}

size_t ElementRepository::Find(const std::string& id, std::vector<ModelElement*>& elements)
{
	remove_copy_if(theElements.begin(), theElements.end(),
		back_inserter(elements), MatchID(id, true));
	return elements.size();
}
