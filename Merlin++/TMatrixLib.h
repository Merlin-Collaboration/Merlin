/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_TMatrixLib
#define _h_TMatrixLib 1

#include "merlin_config.h"
#include <cstdlib>
#include <utility>
#include <cassert>
#include <cmath>
#include <valarray>
#include <iostream>

#include "tblas.h"
#include "utils.h"

namespace TLAS
{

// Defined classes.
template<class T> class Vector;
template<class T> class Matrix;
template<class T> class SubVector;
template<class T> class SubMatrix;
template<class T> class ConstSubVector;
template<class T> class ConstSubMatrix;
struct IdentityMatrix;

// exceptions
class DimensionError
{
};
class RangeError
{
};

/**
 * \class Subscript
 *
 * Single integer subscript to an array object - size_t is unsigned.
 */
typedef size_t Subscript;
typedef size_t Dimension;

/**
 * \class Range
 *
 * A contiguous inclusive range of Subscript values - also unsigned.
 */
typedef std::pair<Subscript, Subscript> Range;

/**
 * Dimensions of a matrix.
 */
struct MatrixDim
{
	Dimension n, m;
	MatrixDim(Dimension nrow, Dimension ncol) :
		n(nrow), m(ncol)
	{
	}
	bool operator==(const MatrixDim& rhs) const
	{
		return n == rhs.n && m == rhs.m;
	}
	bool operator!=(const MatrixDim& rhs) const
	{
		return n != rhs.n || m != rhs.m;
	}

};

inline Dimension length(const Range& r)
{
	return r.second - r.first + 1;
}

inline MatrixDim swap(const MatrixDim& r)
{
	return MatrixDim(r.m, r.n);
}

/**
 * Helper class for constructing identity matrix
 */
struct IdentityMatrix
{
	Dimension dim;
	IdentityMatrix(Dimension n) :
		dim(n)
	{
	}

};

#if !defined(TLAS_NO_RANGE_CHECKING)
inline void equal_dimension(const MatrixDim& d1, const MatrixDim& d2)
{
	if(d1 != d2)
	{
		throw DimensionError();
	}
}
inline void equal_length(Dimension l1, Dimension l2)
{
	if(l1 != l2)
	{
		throw DimensionError();
	}
}
inline void in_range(Subscript i, Dimension l)
{
	if(i >= l)
	{
		throw RangeError();
	}
}
inline void in_range(Range r, Dimension n)
{
	if((r.first >= n) || (r.second >= n))
	{
		throw RangeError();
	}
}

//Pointless function
//inline void valid_dimension(Dimension d1)
//{
//    if(d1<0)
//        throw DimensionError();
//}
#define valid_dimension(dummy) ((void) 0)

#else
#define equal_dimension(dummy1, dummy2) ((void) 0)
#define equal_length(dummy1, dummy2) ((void) 0)
#define in_range(dummy1, dummy2) ((void) 0)
#define valid_dimension(dummy) ((void) 0)
#endif // TLAS_NO_RANGE_CHECKING

// The following typedef is used to put the typedefs into
// the scope of each class requiring them. Use of a base class
// (TMTRX_BASE) which is the correct way to do it is prohibited
// because of a bug in GCC 3.1

#define VEC_TYPE_DEFS \
	typedef std::valarray<T> array_type; \
	typedef std::slice slice_type; \
	typedef std::gslice gslice_type; \
	typedef std::valarray<size_t> index_array_type; \
	typedef T* iterator; \
	typedef const T* const_iterator

/**
 * TMTRX_BASE class,
 */
template<class T>
class TMTRX_BASE
{
public:
	// public typedefs
	typedef T value_type;
	VEC_TYPE_DEFS;
};

template<class T>
class Vector /*: public TMTRX_BASE<T>*/
{

public:

	typedef T value_type;
	VEC_TYPE_DEFS;

	// construction/destruction
	Vector() :
		array()
	{
	}
	explicit Vector(Dimension n) :
		array(n)
	{
	}
	Vector(const T& x, Dimension n) :
		array(x, n)
	{
	}
	Vector(const T* p, Dimension n) :
		array(p, n)
	{
	}
	Vector(const Vector<T>& rhs) :
		array(rhs.array)
	{
	}
	~Vector()
	{
	}

	// copy assignment
	// Non-resize copy (fast)
	const Vector<T>& operator=(const Vector& rhs)
	{
		if(this != &rhs)
		{
			equal_length(array.size(), rhs.array.size());
			array = rhs.array;
		}
		return *this;
	}

	// Re-sizing copy
	const Vector<T>& copy(const Vector& rhs)
	{
		if(this != &rhs)
		{
			if(rhs.array.size() != array.size())
			{
				array.resize(rhs.array.size());
			}
			array = rhs.array;
		}
		return *this;
	}

	const Vector<T>& operator=(const T& val)
	{
		array = val;
		return *this;
	}

	// vector size.
	Dimension size() const
	{
		return array.size();
	}

	void redim(Dimension d)
	{
		array.resize(d);
	}

	// element access
	// We provide both C-like [] and FORTRAN-like ()
	// indexing. Indices start at 0 (C-like).
	T& operator[](Subscript i)
	{
		in_range(i, array.size());
		return array[i];
	}

	T& operator()(Subscript i)
	{
		in_range(i, array.size());
		return array[i];
	}

	T operator[](Subscript i) const
	{
		in_range(i, size());
		return array[i];
	}

	T operator()(Subscript i) const
	{
		in_range(i, size());
		return array[i];
	}

	SubVector<T> operator()(const Range& r)
	{
		in_range(r, size());
		return SubVector<T>(array, slice_type(r.first, length(r), 1));
	}

	SubVector<T> operator[](const Range& r)
	{
		in_range(r, size());
		return SubVector<T>(array, slice_type(r.first, length(r), 1));
	}

	ConstSubVector<T> operator()(const Range& r) const
	{
		in_range(r, size());
		return ConstSubVector<T>(array, slice_type(r.first, length(r), 1));
	}

	ConstSubVector<T> operator[](const Range& r) const
	{
		in_range(r, size());
		return ConstSubVector<T>(array, slice_type(r.first, length(r), 1));
	}

	// logical comparisons
	bool operator==(const Vector<T>& v) const
	{
		equal_length(array.size(), v.array.size());
		for(size_t i = 0; i < array.size(); i++)
			if(array[i] != v.array[i])
			//if(!fequal(array[i],v.array[i]))
			{
				return false;
			}
		return true;
	}

	bool operator!=(const Vector<T>& v) const
	{
		return !operator==(v);
	}

	// s operations
	const Vector<T>& operator+=(const T& s)
	{
		array += s;
		return *this;
	}
	const Vector<T>& operator-=(const T& s)
	{
		array -= s;
		return *this;
	}
	const Vector<T>& operator*=(const T& s)
	{
		array *= s;
		return *this;
	}
	const Vector<T>& operator/=(const T& s)
	{
		array /= s;
		return *this;
	}
       

	// vector operations
	const Vector<T>& operator+=(const Vector<T>& u)
	{
		equal_length(size(), u.size());
		array += u.array;
		return *this;
	}
	const Vector<T>& operator-=(const Vector<T>& u)
	{
		equal_length(size(), u.size());
		array -= u.array;
		return *this;
	}
	const Vector<T>& operator*=(const Vector<T>& u)
	{
		equal_length(size(), u.size());
		array *= u.array;
		return *this;
	}
	const Vector<T>& operator/=(const Vector<T>& u)
	{
		equal_length(size(), u.size());
		array /= u.array;
		return *this;
	}

	// iterator constructors

	iterator begin()
	{
		return std::begin(array);
	}

	iterator end()
	{
		return std::end(array);
	}

	const_iterator begin() const
	{
		return const_cast<Vector<T>*>(this)->begin();
	}

	const_iterator end() const
	{
		return const_cast<Vector<T>*>(this)->end();
	}

private:

	// private constructor for SubVector
	Vector(const array_type& sla) :
		array(sla)
	{
	}

	array_type array;

	friend class SubVector<T>;
	friend class ConstSubVector<T>;
	friend class Matrix<T>;
};

// binary operations for Vector
template<class T>
Vector<T> operator+(const Vector<T>& u, const Vector<T>& v)
{
	Vector<T> rv(u);
	return rv += v;
}

template<class T>
Vector<T> operator-(const Vector<T>& u, const Vector<T>& v)
{
	Vector<T> rv(u);
	return rv -= v;
}

template<class T>
Vector<T> operator/(const Vector<T>& u, const Vector<T>& v)
{
	Vector<T> rv(u);
	return rv /= v;
}

/**
 * Dot (inner) product
 */
template<class T>
T operator*(const Vector<T>& u, const Vector<T>& v)
{
	return tblas1::tdot(u, v, T(0));
}
template<class T>
T operator*(const SubVector<T>& u, const SubVector<T>& v) // RJB added
{
	return tblas1::tdot(u, v, T(0));
}

// scalar-vector algebra

template<class T>
Vector<T> operator+(T s, const Vector<T>& v)
{
	Vector<T> rv(v);
	return rv += s;
}

template<class T>
Vector<T> operator+(const Vector<T>& v, T s)
{
	Vector<T> rv(v);
	return rv += s;
}

template<class T>
Vector<T> operator-(T s, const Vector<T>& v)
{
	Vector<T> rv(v);
	rv *= -1;
	return rv += s;
}

template<class T>
Vector<T> operator-(const Vector<T>& v, T s)
{
	Vector<T> rv(v);
	return rv += s;
}

template<class T>
Vector<T> operator*(T s, const Vector<T>& v)
{
	Vector<T> rv(v);
	return rv *= s;
}

template<class T>
Vector<T> operator*(const Vector<T>& v, T s)
{
	Vector<T> rv(v);
	return rv *= s;
}

template<class T>
Vector<T> operator/(const Vector<T>& v, T s)
{
	Vector<T> rv(v);
	return rv /= s;
}

/*
 * template<class T> SubVector
 * A single (1-dimensional) slice_type (vector)
 * of a Matrix or Vector. SubVector objects have
 * reference semantics to their associated
 * matrix. They cannot be constructed or assigned
 * (except by a Matrix). They are intended as a
 * helper class for accessing rows and columns of
 * Matrices.
 */

template<class T>
class ConstSubVector /* : public TMTRX_BASE<T> */
{
public:

	VEC_TYPE_DEFS;

	// public typedefs
	typedef T value_type;
	typedef Vector<T> vector_type;
	typedef Matrix<T> matrix_type;

	// vector subscripting.
	// () and [] are supported.

	T operator()(Subscript i) const
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	T operator[](Subscript i) const
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	/**
	 * Conversion to a Vector
	 */
	operator Vector<T> () const
	{
		return Vector<T>(array[sl]);
	}

	/**
	 * Length of the vector
	 */
	Dimension size() const
	{
		return sl.size();
	}

	// copying  allowed
	ConstSubVector(const ConstSubVector<T>& v) :
		sl(v.sl), array(v.array)
	{
	}
	ConstSubVector<T>& operator=(const ConstSubVector<T>&)
	{
		return *this;
	}
private:
	// Constructor only called by Matrix<T> (via friendship)
	ConstSubVector(const array_type& ar, const slice_type& s) :
		sl(s), array(ar)
	{
	}

	slice_type sl;
	const array_type& array;

	friend class Matrix<T>;
	friend class Vector<T>;
	friend class SubVector<T>;
};

template<class T>
class SubVector /*: public TMTRX_BASE<T> */
{
public:

	typedef T value_type;
	VEC_TYPE_DEFS;

	// public typedefs
	typedef Vector<T> vector_type;
	typedef Matrix<T> matrix_type;

        T sumsquared() // added RJB
        {
               T s=0;
               for(size_t i=0;i<sl.size();i++) s += pow(array[sl.start()+i*sl.stride()] ,2);  
               return s;
        }
        
        T sumfabs() // added RJB
        {
               T s=0;
               for(size_t i=0;i<sl.size();i++) s += fabs(array[sl.start()+i*sl.stride()] );  
               return s;
        }

        size_t locmaxabs() // added RJB
        {
               size_t imax=0;

               double amax = fabs(array[sl.start()]);
               if(sl.size()>1){
                     for (size_t i=1; i < sl.size();i++) {
                            if(fabs(array[sl.start()+i*sl.stride()]) > amax) {
                            amax = fabs(array[sl.start() + i * sl.stride()]);
                            imax = i;
                        }
                      }
               }

         return imax;
         }

	// vector subscripting.
	// () and [] are supported.

	T& operator()(Subscript i)
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	const T& operator()(Subscript i) const
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	T& operator[](Subscript i)
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	const T& operator[](Subscript i) const
	{
		in_range(i, sl.size());
		return array[sl.start() + i * sl.stride()];
	}

	/**
	 * Assignment to a Vector.
	 */
	SubVector<T>& operator=(const Vector<T>& v)
	{
		equal_length(sl.size(), v.size());
		array[sl] = v.array;
		return *this;
	}

	// assign all elements to a constant scaler
	SubVector<T>& operator=(const T& s)
	{
		array[sl] = s;
		return *this;
	}

	SubVector<T>& operator=(const ConstSubVector<T>& sv)
	{
		equal_length(sl.size(), sv.sl.size());
		array[sl] = sv.array[sv.sl];
		return *this;
	}

	// arithmetic assignment operations (scaler)
	SubVector<T>& operator+=(const T& s)
	{
		array[sl] += array_type(s, sl.size());
		return *this;
	}
	SubVector<T>& operator-=(const T& s)
	{
		array[sl] -= array_type(s, sl.size());
		return *this;
	}
	SubVector<T>& operator*=(const T& s)
	{
		array[sl] *= array_type(s, sl.size());
		return *this;
	}
        Vector<T> operator*(const T& s) const // RJB added
        {
                Vector<T> v(*this);
                v *= s;
                return v;
        }

	SubVector<T>& operator/=(const T& s)
	{
		array[sl] /= array_type(s, sl.size());
		return *this;
	}

	SubVector<T>& operator+=(const Vector<T>& v)
	{
		equal_length(sl.size(), v.size());
		array[sl] += v.array;
		return *this;
	}
	SubVector<T>& operator-=(const Vector<T>& v)
	{
		equal_length(sl.size(), v.size());
		array[sl] -= v.array;
		return *this;
	}
	SubVector<T>& operator*=(const Vector<T>& v)
	{
		equal_length(sl.size(), v.size());
		array[sl] *= v.array;
		return *this;
	}
	SubVector<T>& operator/=(const Vector<T>& v)
	{
		equal_length(sl.size(), v.size());
		array[sl] /= v.array;
		return *this;
	}
	Vector<T> operator/(const T& s) // added RJB
        {   
        Vector<T> rv(*this);
	return rv /s;
        }

	/**
	 * Conversion to a Vector
	 */
	operator Vector<T> () const
	{
		return Vector<T>(array[sl]);
	}

	/**
	 * Convert to a ConstSubVector
	 */
	operator ConstSubVector<T> () const
	{
		return ConstSubVector<T>(array, sl);
	}

	/**
	 * Length of the vector
	 */
	Dimension size() const
	{
		return sl.size();
	}

private:

	/**
	 * Constructor only called by Matrix<T> (via friendship)
	 */
	SubVector(array_type& ar, const slice_type& s) :
		sl(s), array(ar)
	{
	}

	/**
	 * Assignment not allowed.
	 */
	SubVector<T>& operator=(const SubVector<T>&)
	{
		return *this;
	}

	slice_type sl;
	array_type& array;

	friend class Matrix<T>;
	friend class Vector<T>;
};
        template<class T> Vector<T> operator*(const T& s, const SubVector<T> & V) // RJB added
        {
            return(V*s);
        }
/**
 * template<class T> SubMatrix
 * A contiguous (2-dimensional) block (sub-matrix)
 * of a Matrix. SubMatrix objects have
 * reference semantics to their associated
 * matrix. They cannot be constructed or assigned
 * (except by a Matrix). They are intended as a
 * helper class for accessing sub-sections of
 * Matrices.
 */

template<class T>
class ConstSubMatrix /* : public TMTRX_BASE<T> */
{

	VEC_TYPE_DEFS;

	/**
	 * Private function to calculate index
	 */
	Subscript index(Subscript i, Subscript j) const
	{
		in_range(i, gl.size()[0]);
		in_range(j, gl.size()[1]);
		return gl.start() + i * gl.stride()[0] + j * gl.stride()[1];
	}

public:
	// public typedefs
	typedef T value_type;
	typedef Vector<T> vector_type;
	typedef Matrix<T> matrix_type;

	/**
	 * Vector subscripting
	 */
	T operator()(Subscript i, Subscript j) const
	{
		return array[index(i, j)];
	}

	/**
	 * Conversion to Matrix
	 */
	operator Matrix<T>() const
	{
		return Matrix<T>(gl.size()[0], gl.size()[1], array[gl]);
	}

	// dimension accessors
	Dimension nrows() const
	{
		return gl.size()[0];
	}
	Dimension ncols() const
	{
		return gl.size()[1];
	}
	MatrixDim dim() const
	{
		return MatrixDim(gl.size()[0], gl.size()[1]);
	}

private:

	/**
	 * Constructor only called by Matrix<T> (via friendship)
	 */
	ConstSubMatrix(const array_type& ar, const gslice_type& s) :
		gl(s), array(ar)
	{
	}

	// copying not allowed
	ConstSubMatrix(const ConstSubMatrix<T>& t) :
		gl(t.gl), array(t.array)
	{
	}                                // needed for explicit instantiation
	ConstSubMatrix<T>& operator=(const ConstSubMatrix<T>&)
	{
		return *this;
	}

	gslice_type gl;
	const array_type& array;

	friend class Matrix<T>;
	friend class SubMatrix<T>;
};

template<class T>
class SubMatrix /* : public TMTRX_BASE<T> */
{

	VEC_TYPE_DEFS;

	/**
	 * Private function to calculate index
	 */
	Subscript index(Subscript i, Subscript j) const
	{
		in_range(i, gl.size()[0]);
		in_range(j, gl.size()[1]);
		return gl.start() + i * gl.stride()[0] + j * gl.stride()[1];
	}

public:
	// public typedefs
	typedef T value_type;
	typedef Vector<T> vector_type;
	typedef Matrix<T> matrix_type;

	// vector subscripting
	T& operator()(Subscript i, Subscript j)
	{
		return array[index(i, j)];
	}

	const T& operator()(Subscript i, Subscript j) const
	{
		return array[index(i, j)];
	}

	/**
	 * Assignment to a Vector.
	 */
	SubMatrix<T>& operator=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array[gl] = m.array;
		return *this;
	}

	// assign all elements to a constant scaler
	SubMatrix<T>& operator=(const T& s)
	{
		array[gl] = s;
		return *this;
	}

	SubMatrix<T>& operator=(const ConstSubMatrix<T>& sm)
	{
		equal_dimension(dim(), sm.dim());
		array[gl] = sm.array[sm.gl];
		return *this;
	}

	// arithmetic assignment operations (scaler)
	SubMatrix<T>& operator+=(const T& s)
	{
		array[gl] += array_type(s, gl.size()[0] * gl.size()[1]);
		return *this;
	}
	SubMatrix<T>& operator-=(const T& s)
	{
		array[gl] -= array_type(s, gl.size()[0] * gl.size()[1]);
		return *this;
	}
	SubMatrix<T>& operator*=(const T& s)
	{
		array[gl] *= array_type(s, gl.size()[0] * gl.size()[1]);
		return *this;
	}
	SubMatrix<T>& operator/=(const T& s)
	{
		array[gl] /= array_type(s, gl.size()[0] * gl.size()[1]);
		return *this;
	}

	// arithmetic assignment operations (matrix)
	SubMatrix<T>& operator+=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array[gl] += m.array;
		return *this;
	}
	SubMatrix<T>& operator-=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array[gl] -= m.array;
		return *this;
	}
	SubMatrix<T>& operator*=(const Matrix<T>& m)
	{
		equal_length(ncols(), m.nrows());
		array[gl] *= m.array;
		return *this;
	}
	SubMatrix<T>& operator/=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array[gl] /= m.array;
		return *this;
	}

	// conversion to Matrix
	operator Matrix<T>() const
	{
		return Matrix<T>(gl.size()[0], gl.size()[1], array[gl]);
	}

	operator ConstSubMatrix<T>() const
	{
		return ConstSubMatrix<T>(array, gl);
	}

	// dimension accessors
	Dimension nrows() const
	{
		return gl.size()[0];
	}
	Dimension ncols() const
	{
		return gl.size()[1];
	}
	MatrixDim dim() const
	{
		return MatrixDim(gl.size()[0], gl.size()[1]);
	}

private:

	/**
	 * Constructor only called by Matrix<T> (via friendship)
	 */
	SubMatrix(array_type& ar, const gslice_type& s) :
		gl(s), array(ar)
	{
	}

	/**
	 * assignment not allowed.
	 */
	SubMatrix<T>& operator=(const SubMatrix<T>&)
	{
		return *this;
	}

	gslice_type gl;
	array_type& array;

	friend class Matrix<T>;
};

template<class T>
class Matrix /* : public TMTRX_BASE<T> */
{
public:
	typedef T value_type;
	typedef std::valarray<T> array_type;
	typedef std::slice slice_type;
	typedef std::gslice gslice_type;
	typedef std::valarray<size_t> index_array_type;

private:
	/**
	 * Private function to calculate index
	 */
	Subscript index(Subscript i, Subscript j) const
	{
		in_range(i, nr);
		in_range(j, nc);
		return i * nc + j;
	}

public:

	// iterators
	typedef T* iterator;
	typedef const T* const_iterator;

	// construction (arb. dimension)
	Matrix() :
		nr(0), nc(0), array(0)
	{
	}

	Matrix(Dimension rows, Dimension cols) :
		nr(rows), nc(cols), array(rows * cols)
	{
		valid_dimension(rows);
		valid_dimension(cols);
	}
	Matrix(Dimension rows, Dimension cols, const T& val) :
		nr(rows), nc(cols), array(val, rows * cols)
	{
		valid_dimension(rows);
		valid_dimension(cols);
	}

	/**
	 * Construction of an identity matrix
	 */
	Matrix(IdentityMatrix I) :
		nr(I.dim), nc(I.dim), array(T(0), I.dim * I.dim)
	{
		valid_dimension(I.dim);
		// sets diagonal to 1
		array[slice_type(0, I.dim, I.dim + 1)] = T(1);
	}

	// copy another matrix type
	// if conversion from U to T allowed
	template<class U>
	Matrix(const Matrix<U>& rhs) :
		nr(rhs.nrows()), nc(rhs.ncols()), array(rhs.nrows() * rhs.ncols())
	{
		iterator q = begin();
		for(typename Matrix<U>::const_iterator rhsp = rhs.begin(); rhsp != rhs.end(); rhsp++, q++)
		{
			*q = T(*rhsp);
		}
		//std::copy(rhs.begin(),rhs.end(),begin());
	}

	// Need an explicit version that compiler can recognize as a copy constructor
	Matrix(const Matrix<T>& rhs) :
		nr(rhs.nrows()), nc(rhs.ncols()), array(rhs.nrows() * rhs.ncols())
	{
		iterator q = begin();
		for(const_iterator rhsp = rhs.begin(); rhsp != rhs.end(); rhsp++, q++)
		{
			*q = T(*rhsp);
		}
	}

	// destruction
	~Matrix()
	{
		/* nothing to do here */
	}

	// copy assignment
	Matrix<T>& operator=(const Matrix<T>& m)
	{
		if(this != &m)
		{
			equal_dimension(dim(), m.dim());
			array = m.array;
		}
		return *this;
	}

	Matrix<T>& operator=(const T& val)
	{
		array = val;
		return *this;
	}

	Matrix<T>& operator=(const ConstSubMatrix<T>& sm)
	{
		assert(&array != &(sm.array)); // that's a no-no
		equal_dimension(dim(), sm.dim());
		array = sm.array[sm.gl];
		return *this;
	}

	Matrix<T>& copy(const Matrix<T>& rhs)
	{
		redim(rhs.dim());
		return operator=(rhs);
	}

	Matrix<T>& copy(const ConstSubMatrix<T>& rhs)
	{
		redim(rhs.dim());
		return operator=(rhs);
	}

	Matrix<T>& copy(IdentityMatrix I)
	{
		redim(I.dim, I.dim);
		array = T(0);
		// sets diagonal to 1
		array[slice_type(0, I.dim, I.dim + 1)] = T(1);
		return *this;
	}

	// dimension accessors
	Dimension nrows() const
	{
		return nr;
	}
	Dimension ncols() const
	{
		return nc;
	}
	MatrixDim dim() const
	{
		return MatrixDim(nr, nc);
	}

	void redim(Dimension n, Dimension m)
	{
		nr = n;
		nc = m;
		array.resize(nr * nc);
	}

	void redim(MatrixDim md)
	{
		nr = md.n;
		nc = md.m;
		array.resize(nr * nc);
	}

	// matrix part accessors
	T& operator()(Subscript i, Subscript j)
	{
		return array[index(i, j)];
	}
	T operator()(Subscript i, Subscript j) const
	{
		return array[index(i, j)];
	}

	/**
	 * Return a ranged row vector
	 * @return Ranged row vector
	 */
	SubVector<T> operator()(Subscript i, Range col_r)
	{
		in_range(i, nr);
		in_range(col_r, nc);
		return SubVector<T>(array, slice_type(nc * i + col_r.first, length(col_r), 1));
	}

	/**
	 * Return a ranged column vector
	 * @return Ranged column vector
	 */
	SubVector<T> operator()(Range row_r, Subscript j)
	{
		in_range(row_r, nr);
		in_range(j, nc);
		return SubVector<T>(array, slice_type(j + row_r.first * nc, length(row_r), nc));
	}

	/**
	 * Return a sub-matrix
	 * @return Sub-matrix
	 */
	SubMatrix<T> operator()(const Range& row_r, const Range& col_r)
	{
		in_range(row_r, nr);
		in_range(col_r, nc);
		index_array_type lengths(2);
		index_array_type strides(2);
		lengths[0] = length(row_r);
		lengths[1] = length(col_r);
		strides[0] = nc;
		strides[1] = 1;
		return SubMatrix<T>(array, gslice_type(index(row_r.first, col_r.first), lengths, strides));
	}

	/**
	 * Return a row.
	 * @return A row
	 */
	SubVector<T> row(Subscript i)
	{
		in_range(i, nr);
		return SubVector<T>(array, slice_type(i * nc, nc, 1));
	}

	/**
	 * Return a column
	 * @return A column
	 */
	SubVector<T> column(Subscript j)
	{
		in_range(j, nc);
		return SubVector<T>(array, slice_type(j, nr, nc));
	}

	/**
	 * Return a ranged row vector
	 * @return Ranged row vector
	 */
	ConstSubVector<T> operator()(Subscript i, Range col_r) const
	{
		in_range(i, nr);
		in_range(col_r, nc);
		return ConstSubVector<T>(array, slice_type(nc * i + col_r.first, length(col_r), 1));
	}

	/**
	 * Return a ranged column vector
	 * @return Ranged column vector
	 */
	ConstSubVector<T> operator()(Range row_r, Subscript j) const
	{
		in_range(row_r, nr);
		in_range(j, nc);
		return ConstSubVector<T>(array, slice_type(j + row_r.first * nc, length(row_r), nc));
	}

	/**
	 * Return a sub-matrix
	 * @return Sub-matrix
	 */
	ConstSubMatrix<T> operator()(const Range& row_r, const Range& col_r) const
	{
		in_range(row_r, nr);
		in_range(col_r, nc);
		index_array_type lengths(2);
		index_array_type strides(2);
		lengths[0] = length(row_r);
		lengths[1] = length(col_r);
		strides[0] = nc;
		strides[1] = 1;
		return ConstSubMatrix<T>(array, gslice_type(index(row_r.first, col_r.first), lengths, strides));
	}

	/**
	 * Return a row.
	 * @return Row
	 */
	ConstSubVector<T> row(Subscript i) const
	{
		in_range(i, nr);
		return ConstSubVector<T>(array, slice_type(i * nc, nc, 1));
	}

	/**
	 * Return a column
	 * @return Column
	 */
	ConstSubVector<T> column(Subscript j) const
	{
		in_range(j, nc);
		return ConstSubVector<T>(array, slice_type(j, nr, nc));
	}

	// C array style accessors
	SubVector<T> operator[](Subscript i)
	{
		return row(i);
	}

	ConstSubVector<T> operator[](Subscript i) const
	{
		return row(i);
	}

	// arithmetic assignment operations (s)
	Matrix<T>& operator+=(const T& val)
	{
		array += val;
		return *this;
	}

	Matrix<T>& operator-=(const T& val)
	{
		array -= val;
		return *this;
	}

	Matrix<T>& operator*=(const T& val)
	{
		array *= val;
		return *this;
	}

	Matrix<T>& operator/=(const T& val)
	{
		array /= val;
		return *this;
	}

	// arithmetic assignment operations (matrix)

	Matrix<T>& operator+=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array += m.array;
		return *this;
	}

	Matrix<T>& operator-=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array -= m.array;
		return *this;
	}

	Matrix<T>& operator*=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array *= m.array;
		return *this;
	}

	Matrix<T>& operator/=(const Matrix<T>& m)
	{
		equal_dimension(dim(), m.dim());
		array /= m.array;
		return *this;
	}

	// stl-like iterator access to raw data
	iterator begin()
	{
		return std::begin(array);
	}
	iterator end()
	{
		return std::end(array);
	}
	const_iterator begin() const
	{
		return const_cast<Matrix<T>*>(this)->begin();
	}
	const_iterator end() const
	{
		return const_cast<Matrix<T>*>(this)->end();
	}

private:

	Dimension nr;
	Dimension nc;
	array_type array;

	friend class SubVector<T>;
	friend class SubMatrix<T>;
	friend class ConstSubMatrix<T>;

	/**
	 * Private constructor for SubMatrix
	 */
	Matrix(Dimension rows, Dimension cols, const array_type& gsa) :
		nr(rows), nc(cols), array(gsa)
	{
	}
};

// binary operations for Matrix<T>
template<class T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T> B)
{
	Matrix<T> rv(A);
	return rv += B;
}
template<class T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T> B)
{
	Matrix<T> rv(A);
	return rv -= B;
}
template<class T>
Matrix<T> operator/(const Matrix<T>& A, const Matrix<T> B)
{
	Matrix<T> rv(A);
	return rv /= B;
}

template<class T>
Matrix<T> operator*(const Matrix<T>& A, T s)
{
	Matrix<T> rv(A);
	return rv *= s;
}

template<class T>
Matrix<T> operator*(T s, const Matrix<T>& A)
{
	Matrix<T> rv(A);
	return rv *= s;
}

/**
 * Matrix multiplication
 */
template<class T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T> B)
{
	equal_dimension(A.dim(), swap(B.dim()));
	Matrix<T> C(A.nrows(), B.ncols(), T(0));
	tblas3::tgemm(false, false, T(1), A, B, T(1), C);
	return C;
}

/**
 * Matrix-vector multiplication
 */
template<class T>
Vector<T> operator*(const Matrix<T>& M, const Vector<T> V)
{
	equal_length(M.ncols(), V.size());
	Vector<T> U(T(0), M.nrows());
	tblas2::tgemv(false, T(1), M, V, T(1), U);
	return U;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& M)
{
	for(size_t i = 0; i < M.ncols(); i++)
	{
		for(size_t j = 0; j < M.nrows(); j++)
		{
			os << M(i, j) << " ";
		}
		os << std::endl;
	}
	return os;
}

extern template class Vector<double>;
extern template class Matrix<double>;
extern template class SubVector<double>;
extern template class SubMatrix<double>;
extern template class ConstSubVector<double>;
extern template class ConstSubMatrix<double>;

} // end of namespace TLAS

#endif  // _h_TMatrixLib
