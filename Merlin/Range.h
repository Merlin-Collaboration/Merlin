/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 * This file is derived from software bearing the copyright notice: Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED.
 */

#ifndef NumericalRange_h
#define NumericalRange_h 1

#include <algorithm>
#include <functional>
#include "utils.h"

class RangeBase
{
public:
	typedef enum
	{
		ok,
		belowlower,
		aboveupper

	} Result;
	virtual ~RangeBase()
	{
	}
};

/**
 *	Represents an allowed contiguous range of a floating
 *	point number.
 */
struct UnboundedRange {};

template<class T, class C = std::less<T> >
class NumericalRange: public RangeBase
{
public:

	/**
	 *	Constructor taking the range (lo,hi>
	 */
	NumericalRange(const T& lo, const T& hi);

	/**
	 *    Construct an unbounded range
	 */
	NumericalRange(UnboundedRange) :
		lower(1), upper(-1)
	{
	}

	/**
	 *    Construct a fixed point (zero range)
	 */
	explicit NumericalRange(const T& fp) :
		lower(fp), upper(fp)
	{
	}

	/**
	 *	Returns true if this range is unbounded (i.e. represents
	 *	+/-infinity).
	 *
	 *	@retval true If range unbounded (i.e. \f$\pm\infty\f$)
	 *	@retval false If range bounded
	 */
	bool IsUnbounded() const;
	bool IsFixedPoint() const
	{
		return lower == upper;
	}

	/**
	 *	Returns true if lower<=x<=upper.
	 *
	 *	@retval true If \f$ \mathrm{lower} \le x \le \mathrm{upper} \f$
	 */
	bool operator ()(const T& x) const;
	bool operator ==(const NumericalRange& rhs) const;
	bool operator !=(const NumericalRange& rhs) const;

	/**
	 *	Checks id x is within the range. Returns belowlower, ok or
	 *	aboveupper.
	 */
	RangeBase::Result Check(const T& x) const;

	/**
	 *	The lowerimum value of the range.
	 */
	T lower;

	/**
	 *	The upperimum value of the range.
	 */
	T upper;
};

template<class T, class C>
inline NumericalRange<T, C>::NumericalRange(const T& lo, const T& hi) :
	lower(lo), upper(hi)
{
	if(C()(upper, lower))
	{
		std::swap(lower, upper);
	}
}

template<class T, class C>
inline bool NumericalRange<T, C>::IsUnbounded() const
{
	return fequal(lower, 1.0) && fequal(upper, -1.0);
}

template<class T, class C>
inline bool NumericalRange<T, C>::operator ()(const T& x) const
{
	if(IsUnbounded())
	{
		return true;
	}
	else
	{
		return x >= lower && x <= upper;
	}
}

template<class T, class C>
inline bool NumericalRange<T, C>::operator ==(const NumericalRange& rhs) const
{
	return lower == rhs.lower && upper == rhs.upper;
}

template<class T, class C>
inline bool NumericalRange<T, C>::operator !=(const NumericalRange& rhs) const
{
	return lower != rhs.lower || upper != rhs.upper;
}

template<class T, class C>
inline RangeBase::Result NumericalRange<T, C>::Check(const T& x) const
{
//dk    if(isUnbounded())
	if(IsUnbounded())
	{
		return ok;
	}
	else if(x < lower)
	{
		return belowlower;
	}
	else if(x <= upper)
	{
		return ok;
	}
	else
	{
		return aboveupper;
	}
}

// float numerical range
typedef NumericalRange<double> FloatRange;

#endif
