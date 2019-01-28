/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _CollimatorAperture_h_
#define _CollimatorAperture_h_

#include "SimpleApertures.h"
#include "MaterialDatabase.h"
#include "Material.h"
#include <iostream>

/**********************************************************************
 *
 *	A collimator jaw, aligned to the beam orbit and beta function changes
 *   This does NOT have jaw flatness errors
 *
 **********************************************************************/

class CollimatorAperture: public RectangularAperture
{
protected:
	double alpha;
	double CollimatorLength;
	double jaw_length;
	double x_offset_entry, y_offset_entry;

//Add jaw parameters at exit as well
	double x_offset_exit, y_offset_exit;
	double w_exit, h_exit;
	double cosalpha;
	double sinalpha;

public:
	CollimatorAperture(double w, double h, double t, Material* m, double length, double x_offset_entry = 0.0, double
		y_offset_entry = 0.0);

	void SetExitWidth(double);  //Horizontal
	void SetExitHeight(double); //Vertical
	void SetExitXOffset(double);    //Horizontal
	void SetExitYOffset(double);    //Vertical

	double GetFullEntranceHeight() const;
	double GetFullEntranceWidth() const;

	double GetFullExitHeight() const;
	double GetFullExitWidth() const;

	double GetEntranceXOffset() const;
	double GetEntranceYOffset() const;

	double GetExitXOffset() const;
	double GetExitYOffset() const;

	double GetCollimatorTilt() const;

	virtual bool PointInside(double x, double y, double z) const;
};

/**********************************************************************
 *
 *	A collimator jaw, unaligned to the beam orbit or beta function changes
 *   This does NOT have jaw flatness errors
 *
 **********************************************************************/

class UnalignedCollimatorAperture: public CollimatorAperture
{
public:
	UnalignedCollimatorAperture(double w, double h, double t, Material* m, double length, double x_offset_entry = 0.0,
		double y_offset_entry = 0.0);

	bool PointInside(double x, double y, double z) const;
};

/**
 *	A collimator jaw, aligned to the beam orbit and beta function changes
 *  This has jaw flatness errors
 */

class CollimatorApertureWithErrors: public CollimatorAperture
{
	double ApertureError;
	bool PointInside(double x, double y, double z) const;
};

/**
 *	A collimator jaw, unaligned to the beam orbit or beta function changes
 *  This has jaw flatness errors
 */

class UnalignedCollimatorApertureWithErrors: public UnalignedCollimatorAperture
{
	double ApertureError;
	bool PointInside(double x, double y, double z) const;
};

/**
 *	A collimator jaw, unaligned to the beam orbit or beta function changes
 *  This does NOT have jaw flatness errors
 */

class OneSidedUnalignedCollimatorAperture: public CollimatorAperture
{
public:
	OneSidedUnalignedCollimatorAperture(double w, double h, double t, Material* m, double length, double
		x_offset_entry = 0.0, double y_offset_entry = 0.0);

	bool PointInside(double x, double y, double z) const;
	bool PositiveSide;
	void SetJawSide(bool);
};

#endif
