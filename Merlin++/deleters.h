/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef deleters_h
#define deleters_h 1

#include <utility>

/**
 *	Function object for deleting containers of pointers.
 */

template<class T>
class deleter
{
public:
	/**
	 * Deletes the pointer p.
	 */
	void operator ()(T* p);
};

/**
 *	Function object for deleting associative containers
 *	whose value types are pointers.
 */

template<class key, class val>
class map_deleter
{
public:
	void operator ()(std::pair<key, val*>& arg);
};

template<class T>
inline void deleter<T>::operator ()(T* p)
{
	if(p)
	{
		delete p;
	}
}

template<class key, class val>
inline void map_deleter<key, val>::operator ()(std::pair<key, val*>& arg)
{
	if(arg.second)
	{
		delete arg.second;
	}
}

#endif
