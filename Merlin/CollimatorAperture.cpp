/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "CollimatorAperture.h"
#include "MADInterface.h"
#include "RandomNG.h"

/**********************************************************************
 *
 *	A collimator jaw, aligned to the beam orbit and beta function changes
 *   This does NOT have jaw flatness errors
 *
 **********************************************************************/
CollimatorAperture::CollimatorAperture(double w, double h, double t, Material* m, double length, double x_off, double
	y_off) :
	RectangularAperture(w, h), alpha(t), CollimatorLength(length), x_offset_entry(x_off), y_offset_entry(y_off),
	cosalpha(cos(-t)), sinalpha(sin(-t))
{
	SetMaterial(m);
	x_offset_exit = 0;
	y_offset_exit = 0;
	w_exit = 0;
	h_exit = 0;
	//cout << cosalpha << "\t" << sinalpha << endl;
}

//Checks if particle is in or outside a defined aperture
inline bool CollimatorAperture::PointInside(double x, double y, double z) const
{
	/*
	   We need to calculate several variables:
	   1: The x,y offsets at the z position of the particle.
	   2: The width and hight of the jaw at the z position of the particle.

	   Lets start with the x and y offsets;
	 */

	//y = mx + c
	//m = dy/dx

	double x_off = (z * (x_offset_entry - x_offset_exit) / CollimatorLength) - x_offset_entry;
	double y_off = (z * (y_offset_entry - y_offset_exit) / CollimatorLength) - y_offset_entry;

	//These will give the jaw width and heights to be used. * 0.5 to convert to half width.
	double x_jaw = (z * (w_exit - GetFullWidth()) / CollimatorLength) + GetFullWidth();
	double y_jaw = (z * (h_exit - GetFullHeight()) / CollimatorLength) + GetFullHeight();

	double x1 = ((x + x_off) * cosalpha) - ((y + y_off) * sinalpha);
	double y1 = ((x + x_off) * sinalpha) + ((y + y_off) * cosalpha);

	return fabs(x1) < (x_jaw / 2) && fabs(y1) < (y_jaw / 2);
}

//Sets the jaw width at the exit of the collimator
void CollimatorAperture::SetExitWidth(double width)
{
	w_exit = width;
}

//Sets the jaw height at the exit of the collimator
void CollimatorAperture::SetExitHeight(double height)
{
	h_exit = height;
}

//Sets the x (horizontal) orbit offset at the exit of the collimator
void CollimatorAperture::SetExitXOffset(double x)
{
	x_offset_exit = x;
}

//Sets the y (vertical) orbit offset at the exit of the collimator
void CollimatorAperture::SetExitYOffset(double y)
{
	y_offset_exit = y;
}

double CollimatorAperture::GetFullEntranceHeight() const
{
	return GetFullHeight();
}

double CollimatorAperture::GetFullEntranceWidth() const
{
	return GetFullWidth();
}

double CollimatorAperture::GetFullExitHeight() const
{
	return h_exit;
}

double CollimatorAperture::GetFullExitWidth() const
{
	return w_exit;
}

double CollimatorAperture::GetEntranceXOffset() const
{
	return x_offset_entry;
}

double CollimatorAperture::GetEntranceYOffset() const
{
	return y_offset_entry;
}

double CollimatorAperture::GetExitXOffset() const
{
	return x_offset_exit;
}

double CollimatorAperture::GetExitYOffset() const
{
	return y_offset_exit;
}

double CollimatorAperture::GetCollimatorTilt() const
{
	return alpha;
}

/**********************************************************************
 *
 *	A collimator jaw, unaligned to the beam orbit or beta function changes
 *   This does NOT have jaw flatness errors
 *	Still aligned to the beam size! Flat jaws, parallel to the beam pipe,
 *	but touching the beam envelope on each side
 *	This is the LHC configuration
 *
 **********************************************************************/
UnalignedCollimatorAperture::UnalignedCollimatorAperture(double w, double h, double t, Material* m, double length,
	double x_off, double y_off) :
	CollimatorAperture(w, h, t, m, length, x_off, y_off)
{
	SetMaterial(m);
}

inline bool UnalignedCollimatorAperture::PointInside(double x, double y, double z) const
{
	double x1 = ((x - x_offset_entry) * cosalpha) - ((y - y_offset_entry) * sinalpha);
	double y1 = ((x - x_offset_entry) * sinalpha) + ((y - y_offset_entry) * cosalpha);

	return fabs(x1) < GetFullWidth() / 2 && fabs(y1) < GetFullHeight() / 2;
}

/**********************************************************************
 *
 *	A collimator jaw, aligned to the beam orbit and beta function changes
 *   This has jaw flatness errors
 *
 **********************************************************************/

inline bool CollimatorApertureWithErrors::PointInside(double x, double y, double z) const
{
	double x_off = (z * (x_offset_entry - x_offset_exit) / CollimatorLength) - x_offset_entry;
	double y_off = (z * (y_offset_entry - y_offset_exit) / CollimatorLength) - y_offset_entry;

	//These will give the jaw width and heights to be used. * 0.5 to convert to half width.
	double x_jaw = (z * (w_exit - GetFullWidth()) / CollimatorLength) + GetFullWidth();
	double y_jaw = (z * (h_exit - GetFullHeight()) / CollimatorLength) + GetFullHeight();

	double x1 = ((x + x_off) * cosalpha) - ((y + y_off) * sinalpha);
	double y1 = ((x + x_off) * sinalpha) + ((y + y_off) * cosalpha);
	return fabs(x1) < (x_jaw / 2) && fabs(y1) < (y_jaw / 2);
	x1 += ApertureError * ((pow(z, 2) / jaw_length) - z);
	y1 += ApertureError * ((pow(z, 2) / jaw_length) - z);

	return fabs(x1) < (x_jaw / 2) && fabs(y1) < (y_jaw / 2);
}

/**********************************************************************
 *
 *	A one sided collimator jaw, unaligned to the beam orbit or beta function changes
 *   This does NOT have jaw flatness errors
 *	Still aligned to the beam size! Flat jaws, parallel to the beam pipe,
 *	but touching the beam envelope on each side
 *	This is the LHC configuration
 *
 **********************************************************************/
OneSidedUnalignedCollimatorAperture::OneSidedUnalignedCollimatorAperture(double w, double h, double t, Material* m,
	double length, double x_off, double y_off) :
	CollimatorAperture(w, h, t, m, length, x_off, y_off), PositiveSide(true)
{
	SetMaterial(m);
}

inline bool OneSidedUnalignedCollimatorAperture::PointInside(double x, double y, double z) const
{
	double x1 = ((x - x_offset_entry) * cosalpha) - ((y - y_offset_entry) * sinalpha);
	double y1 = ((x - x_offset_entry) * sinalpha) + ((y - y_offset_entry) * cosalpha);

	if(PositiveSide)
	{
		return x1 < GetFullWidth() / 2 && fabs(y1) < GetFullHeight() / 2;
	}
	else
	{
		return (-x1) < GetFullWidth() / 2 && fabs(y1) < GetFullHeight() / 2;
	}
}

void OneSidedUnalignedCollimatorAperture::SetJawSide(bool side)
{
	PositiveSide = side;
}
