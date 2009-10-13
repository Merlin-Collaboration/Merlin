/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file NumericalUtils\Range.h
 * last modified 04/26/00  11:36:52
 *
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED. 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */

#ifndef NumericalRange_h
#define NumericalRange_h 1

#include <algorithm>
#include <functional>


class RangeBase
{
public:
    typedef enum {ok,belowlower,aboveupper} Result;
};

//	Represents an allowed contiguous range of a floating
//	point number.
struct UnboundedRange {};

template <class T, class C = std::less<T> >
class NumericalRange : public RangeBase {
public:

    //	Constructor taking the range (lo,hi>
    NumericalRange (const T& lo, const T& hi);

    //    Construct an unbounded range
    NumericalRange (UnboundedRange) : lower(1),upper(-1) {}

    //    Construct a fixed point (zero range)
    explicit NumericalRange (const T& fp) : lower(fp),upper(fp) {}

    //	Returns true if this range is unbouned (i.e. represents
    //	+/-infinity).
    bool IsUnbounded () const;
    bool IsFixedPoint () const { return lower==upper; }

    //	Returns true if lower<=x<=upper.
    bool operator () (const T& x) const;
    bool operator == (const NumericalRange& rhs) const;
    bool operator != (const NumericalRange& rhs) const;

    //	Checks id x is within the range. Returns belowlower, ok or
    //	aboveupper.
    RangeBase::Result Check (const T& x) const;

    //	The lowerimum value of the range.
    T lower;

    //	The upperimum value of the range.
    T upper;
};


template <class T, class C>
inline NumericalRange<T,C>::NumericalRange (const T& lo, const T& hi)
        : lower(lo),upper(hi)
{
    if(C()(upper,lower))
        std::swap(lower,upper);
}

template <class T, class C>
inline bool NumericalRange<T,C>::IsUnbounded () const
{
    return lower==1 && upper==-1;
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator () (const T& x) const
{
    if(IsUnbounded())
        return true;
    else
        return x>=lower && x<=upper;
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator == (const NumericalRange& rhs) const
{
    return lower==rhs.lower && upper==rhs.upper;
}

template <class T, class C>
inline bool NumericalRange<T,C>::operator != (const NumericalRange& rhs) const
{
    return lower!=rhs.lower || upper!=rhs.upper;
}

template <class T, class C>
inline RangeBase::Result NumericalRange<T,C>::Check (const T& x) const
{
//dk    if(isUnbounded())
    if(IsUnbounded())
        return ok;
    else if(x<lower)
        return belowlower;
    else if(x<=upper)
        return ok;
    else
        return aboveupper;
}

// float numerical range
typedef NumericalRange< double > FloatRange;


#endif
