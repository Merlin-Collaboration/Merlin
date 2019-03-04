/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BunchConstructor_h
#define BunchConstructor_h 1

#include "merlin_config.h"

class Bunch;

/**
 *	Abstract factory for constructing a Bunch.
 */

class BunchConstructor
{
public:

	virtual ~BunchConstructor();

	/**
	 *	Constructs a (new) bunch in memory. The bunch index is
	 *	supplied for implementations that support multiple
	 *	bunches (i.e. bunch trains).
	 */
	virtual Bunch* ConstructBunch(int bunchIndex = 0) const = 0;
};

/**
 *	Template class to generate a BunchConstructor which
 *	constructs a Bunch of type B. The Construct() method
 *	always returns a copy of the same (stored) bunch. Class
 *	B must provide a copy constructor.
 */

template<class B>
class StaticBunchCtor: public BunchConstructor
{
public:
	explicit StaticBunchCtor(B* source, bool del = false);

	~StaticBunchCtor();

	/**
	 *	Constructs and returns a copy of the source bunch.
	 */
	virtual Bunch* ConstructBunch(int bunchIndex = 0) const;

	/**
	 *	Sets the source bunch. Set del to true if the bunch is
	 *	to be deleted when the destructor is called.
	 */
	void SetSourceBunch(B* bunch0, bool del = false);

private:

	B* sourceBunch;

	/**
	 *	Set true if the ctor owns the source bunch.
	 */
	bool owns;

	//Copy protection for bunches
	StaticBunchCtor(const StaticBunchCtor& rhs);
	StaticBunchCtor& operator=(const StaticBunchCtor& rhs);
};

inline BunchConstructor::~BunchConstructor()
{
}

template<class B>
inline StaticBunchCtor<B>::StaticBunchCtor(B* source, bool del) :
	sourceBunch(source), owns(del)
{
}

template<class B>
StaticBunchCtor<B>::~StaticBunchCtor()
{
	if(owns && sourceBunch != nullptr)
	{
		delete sourceBunch;
	}
}

template<class B>
Bunch* StaticBunchCtor<B>::ConstructBunch(int bunchIndex) const
{
	return new B(*sourceBunch);
}

template<class B>
void StaticBunchCtor<B>::SetSourceBunch(B* bunch0, bool del)
{
	if(owns && sourceBunch != 0)
	{
		delete sourceBunch;
	}
	sourceBunch = bunch0;
	owns = del;
}

template<class B>
inline StaticBunchCtor<B>* MakeBunchCtor(B* bunch0, bool del = false)
{
	return new StaticBunchCtor<B>(bunch0, del);
}

#endif
