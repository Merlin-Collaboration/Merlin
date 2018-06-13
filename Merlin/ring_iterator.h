/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ring_iterator_h
#define ring_iterator_h 1

#include <iterator>

/**
 * template class ring_iterator
 * An iterator which cyclically iterates through a sequence (container).
 */

template<class C, class I>
class ring_iterator
{
private:

	I curr;   /// current iterator position
	C* cont;  /// pointer to container

public:

	typedef typename C::value_type& reference_type;
	typedef typename C::value_type* pointer_type;

	ring_iterator(C& c, I itor) :
		curr(itor), cont(&c)
	{
	}
	ring_iterator()
	{
	}

	ring_iterator& operator++()
	{
		if((++curr) == cont->end())
		{
			curr = cont->begin();
		}
		return *this;
	}

	ring_iterator operator++(int)
	{
		ring_iterator tmp(*this);
		if((++curr) == cont->end())
		{
			curr = cont->begin();
		}
		return tmp;
	}

	ring_iterator& operator--()
	{
		if(curr == cont->begin())
		{
			curr = cont->end();
		}
		--curr;
		return *this;
	}

	ring_iterator operator--(int)
	{
		ring_iterator tmp(*this);
		if(curr == cont->begin())
		{
			curr = cont->end();
		}
		--curr;
		return tmp;
	}

	bool operator==(const ring_iterator<C, I>& rhs) const
	{
		return curr == rhs.curr;
	}

	bool operator!=(const ring_iterator<C, I>& rhs) const
	{
		return curr != rhs.curr;
	}

	bool operator==(const I& rhs) const
	{
		return curr == rhs;
	}

	bool operator!=(const I& rhs) const
	{
		return curr != rhs;
	}

	reference_type operator*()
	{
		return *curr;
	}

	pointer_type operator->()
	{
		return &*curr;
	}
};

template<class C, class I>
inline ring_iterator<C, I> make_ring_iterator(C& cont, I& iter)
{
	return ring_iterator<C, I>(cont, iter);
}

template<class C>
inline ring_iterator<C, typename C::iterator> make_ring_iterator(C& cont)
{
	return ring_iterator<C, typename C::iterator>(cont, cont.begin());
}

#endif
