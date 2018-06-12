/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MultipoleField_h
#define MultipoleField_h 1

#include "merlin_config.h"
#include <vector>
#include <utility>
#include <iostream>
#include <complex>
typedef std::complex<double> Complex;
#include "Space3D.h"
#include "EMField.h"
#include "utils.h"

/**
 *	A representation of a multipole expansion of a magnetic
 *	field having only a non-zero component of the magnetic
 *	vector potential in the z-direction. The general
 *	expression for the series is give as
 */

//	                  n
//	     By+i Bx =  Sum  B0 (bn-i an) (z/r0)
//	                 i=0

/**
 *	\f[
 *        B_y + i B_x = \sum_{j=0}^n B_0 (b_n - i a_n) (\frac{z}{r_0}
 *	\f]
 *
 *	where z=(x+iy). B0 is the overall scale (in Tesla), and
 *	r0 is an arbitrary defining radius. For Multipole
 *	expansions representing ideal single multipoles (e.g.
 *	dipole, quadrupole etc), then B0 is set to the field
 *	strength at r0, and consequently |bn+ i an|=1.
 */

class MultipoleField: public EMField
{
public:

	typedef std::vector<Complex> TermExpansion;
	typedef TermExpansion::iterator iterator;
	typedef TermExpansion::const_iterator const_iterator;

	/**
	 *	Constructor forming a null (zero) field, but with a
	 *	specified field scale factor (in Tesla).
	 */
	MultipoleField(double scale, size_t nt = 0);

	/**
	 *	Constructor forming a single multipole. The field
	 *	components bn and an are in Tesla , and are normalised
	 *	to the radius r`
	 */
	MultipoleField(int np, double bn, double an, double r0);

	/**
	 *	Constructor forming a single multipole. The field
	 *	component bn is in Tesla , and are normalised to the
	 *	radius r. If skew is true, creates a skew-field`
	 */
	MultipoleField(int np, double bn, double r0, bool skew = false);

	/**
	 *	Constructor forming a single multipole. The field
	 *	component strength bn is in Tesla/meter^n , where n is
	 *	the pole number). If skew is true, creates a skew-field.
	 */
	MultipoleField(int np, double bn, bool skew = false);

	/**
	 *	Returns the current normalising scale for the field (in
	 *	Tesla)
	 *	@return Current normalising scale
	 */
	double GetFieldScale() const;

	/**
	 *	Set the (normalising) scale of the field (units are in
	 *	Tesla).
	 */
	void SetFieldScale(double scale);

	/**
	 *	Returns true if this field is a null field.
	 *	@retval true If the field is a null field
	 *	@retval false If the field is not a null field
	 */
	bool IsNullField() const;

	/**
	 *	Returns the normalised field component Kn for the
	 *	specified magnetic rigidity brho. Kn is defined as
	 *	Dn[By+iBx]/brho, where Dn is the n-th complex derivative
	 *	wrt x+i y, evaluated at (x+iy)=0. If brho is in
	 *	Tesla.meter, Kn has the units 1/meter^(n+1).
	 */
	Complex GetKn(int np, double rigidity) const;

	/**
	 *	Returns the 2-D magnetic field in Tesla as the complex
	 *	number By+i Bx at the point (x,y) (in meters). The
	 *	optional parameter exclude can be used to exclude terms
	 *	lower than or equal to exclude Thus a value of exclude=1
	 *	gives only that field due to the non-linear terms
	 *	(sextupole and above). A Value of -1 (default) includes
	 *	all multipole terms.
	 */
	Complex GetField2D(double x, double y, int exclude = -1) const;

	/**
	 *	Returns the magnetic field at the point x. The time
	 *	variable t is ignored.
	 *	@return Magnetic field at point x
	 */
	virtual Vector3D GetBFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Returns a null vector (pure magnetic field).
	 *	@return Null vector
	 */
	virtual Vector3D GetEFieldAt(const Point3D& x, double t = 0) const;

	/**
	 *	Return the force due to the magnetic field on a particle
	 *	of charge q, with position x and velocity v. The time
	 *	variable t is ignored (static magnetic field).
	 *
	 *	@return Lorentz force
	 */
	virtual Vector3D GetForceAt(const Point3D& x, const Vector3D& v, double q, double t = 0) const;

	/**
	 *	Rotates this field 180 degrees about the local vertical
	 *	access.
	 */
	void RotateY180();

	/**
	 *	Prints an ASCII representation (table) of the field to os.
	 */
	void PrintField(std::ostream& os) const;

	/**
	 *	Add the specified component to the field. The field
	 *	components bn and an are in Tesla , and are normalised
	 *	to the radius r (default = 1 meter).
	 */
	void SetComponent(size_t np, double bn, double an = 0, double r0 = 1);

	/**
	 *	Returns the field (in Tesla) at the radius r0 due to the
	 *	np-th multipole component
	 *	@return Magnetic field strength (T) due to np-th multipole
	 */
	Complex GetComponent(size_t np, double r0 = 1.0) const;

	/**
	 *	Returns the unitless complex coefficient for the np-th
	 *	term (bn+i*an).The coefficient is relative to the
	 *	specified pole radius r0 (default = 1meter).
	 *	@return Complex coefficient for np-th term
	 */
	Complex GetCoefficient(size_t np, double r0 = 1.0) const;

	/**
	 *	Sets the unitless complex coefficient for the np-th term
	 *	(bn+i*an).The coefficient is relative to the specified
	 *	pole radius r0 (default = 1meter).
	 */
	void SetCoefficient(size_t np, const Complex& b, double r0 = 1.0);

	/**
	 * Return the highest non-zero multipole index.
	 * @return Highest non-zero multipole index
	 */
	int HighestMultipole() const;

	/**
	 * Return the lowest non-zero multipole index.
	 * @return Lowest non-zero multipole index
	 */
	int LowestMultipole() const;

private:

	/**
	 *	The scale of the field in Tesla. The coefficients of the
	 *	expansion terms are normalised so that they have no
	 *	field unit.
	 */
	double B0;

	TermExpansion expansion;
};

inline MultipoleField::MultipoleField(double scale, size_t nt) :
	B0(scale), expansion(nt, Complex(0, 0))
{
}

inline bool MultipoleField::IsNullField() const
{
	return fequal(B0, 0.0);
}

inline int MultipoleField::HighestMultipole() const
{
	return expansion.size() - 1;
}

inline int MultipoleField::LowestMultipole() const
{
	size_t n = 0;
	while(expansion[n++] == Complex(0.0))
		;
	return n - 1;
}

#endif
