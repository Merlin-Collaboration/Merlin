/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "CollimatorAperture.h"
#include "MADInterface.h"
#include "RandomNG.h"

CollimatorAperture::CollimatorAperture(double w, double h, double t, double length, double x_off, double y_off) :
	alpha(t), CollimatorLength(length), x_offset_entry(x_off), y_offset_entry(y_off), x_offset_exit(0), y_offset_exit(
		0), w_entrance(w), h_entrance(h), w_exit(0), h_exit(0), cosalpha(cos(-t)), sinalpha(sin(-t))
{
	setRectHalfWidth(w / 2);
	setRectHalfHeight(h / 2);
	setMinDim(min((w / 2), (h / 2)));
	setApertureType("COLLIMATOR");
}

inline bool CollimatorAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	if(x_offset_entry == 0 && y_offset_entry == 0 && x_offset_exit == 0 && y_offset_exit == 0)
	{
		double ax = fabs(x);
		double ay = fabs(y);
		if(ax + ay < minDim)
			return true;
	}
	double x_off = (z * (x_offset_entry - x_offset_exit) / CollimatorLength) - x_offset_entry;
	double y_off = (z * (y_offset_entry - y_offset_exit) / CollimatorLength) - y_offset_entry;

	double x_jaw = (z * (w_exit - GetFullEntranceWidth()) / CollimatorLength) + GetFullEntranceWidth();
	double y_jaw = (z * (h_exit - GetFullEntranceHeight()) / CollimatorLength) + GetFullEntranceHeight();

	double x1 = ((x + x_off) * cosalpha) - ((y + y_off) * sinalpha);
	double y1 = ((x + x_off) * sinalpha) + ((y + y_off) * cosalpha);

	return fabs(x1) * 2 < x_jaw && fabs(y1) * 2 < y_jaw;
}

void CollimatorAperture::SetEntranceWidth(double width)
{
	w_entrance = width;
}

void CollimatorAperture::SetEntranceHeight(double height)
{
	h_entrance = height;
}

void CollimatorAperture::SetExitWidth(double width)
{
	w_exit = width;
}

void CollimatorAperture::SetExitHeight(double height)
{
	h_exit = height;
}

void CollimatorAperture::SetExitXOffset(double x)
{
	x_offset_exit = x;
}

void CollimatorAperture::SetExitYOffset(double y)
{
	y_offset_exit = y;
}

double CollimatorAperture::GetFullEntranceHeight() const
{
	return rectHalfHeight * 2;
}

double CollimatorAperture::GetFullEntranceWidth() const
{
	return rectHalfWidth * 2;
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

double CollimatorAperture::GetCollimatorLength() const
{
	return CollimatorLength;
}

UnalignedCollimatorAperture::UnalignedCollimatorAperture(double w, double h, double t, double length, double x_off,
	double y_off) :
	CollimatorAperture(w, h, t, length, x_off, y_off)
{
	setRectHalfWidth(w / 2);
	setRectHalfHeight(h / 2);
	setMinDim(min((w / 2), (h / 2)));
	setApertureType("COLLIMATOR");
}

inline bool UnalignedCollimatorAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double x1 = ((x - x_offset_entry) * cosalpha) - ((y - y_offset_entry) * sinalpha);
	double y1 = ((x - x_offset_entry) * sinalpha) + ((y - y_offset_entry) * cosalpha);

	return fabs(x1) * 2 < GetFullEntranceWidth() && fabs(y1) * 2 < GetFullEntranceHeight();
}

inline bool CollimatorApertureWithErrors::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double x_off = (z * (x_offset_entry - x_offset_exit) / CollimatorLength) - x_offset_entry;
	double y_off = (z * (y_offset_entry - y_offset_exit) / CollimatorLength) - y_offset_entry;

	double x_jaw = (z * (w_exit - GetFullEntranceWidth()) / CollimatorLength) + GetFullEntranceWidth();
	double y_jaw = (z * (h_exit - GetFullEntranceHeight()) / CollimatorLength) + GetFullEntranceHeight();

	double x1 = ((x + x_off) * cosalpha) - ((y + y_off) * sinalpha);
	double y1 = ((x + x_off) * sinalpha) + ((y + y_off) * cosalpha);
	return fabs(x1) * 2 < x_jaw && fabs(y1) * 2 < y_jaw;
	x1 += ApertureError * ((pow(z, 2) / CollimatorLength) - z);
	y1 += ApertureError * ((pow(z, 2) / CollimatorLength) - z);

	return fabs(x1) * 2 < x_jaw && fabs(y1) * 2 < y_jaw;
}

OneSidedUnalignedCollimatorAperture::OneSidedUnalignedCollimatorAperture(double w, double h, double t, double length,
	double x_off, double y_off, bool side) :
	CollimatorAperture(w, h, t, length, x_off, y_off), JawSide(side)
{
	setRectHalfWidth(w / 2);
	setRectHalfHeight(h / 2);
	setMinDim(min((w / 2), (h / 2)));
	setApertureType("ONE-SIDED COLLIMATOR");
}

inline bool OneSidedUnalignedCollimatorAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	double x1 = ((x - x_offset_entry) * cosalpha) - ((y - y_offset_entry) * sinalpha);
	double y1 = ((x - x_offset_entry) * sinalpha) + ((y - y_offset_entry) * cosalpha);

	if(JawSide)
	{
		return x1 * 2 < GetFullEntranceWidth() && fabs(y1) * 2 < GetFullEntranceHeight();
	}
	else
	{
		return -(x1 * 2) < GetFullEntranceWidth() && fabs(y1) * 2 < GetFullEntranceHeight();
	}
}

void OneSidedUnalignedCollimatorAperture::SetJawSide(bool side)
{
	JawSide = side;
}

bool OneSidedUnalignedCollimatorAperture::GetJawSide()
{
	return JawSide;
}
