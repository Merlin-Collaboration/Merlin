/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _CollimatorAperture_h_
#define _CollimatorAperture_h_

#include "Aperture.h"
#include "MaterialDatabase.h"
#include "Material.h"
#include <iostream>

/**
 * Collimator aperture class
 */
class CollimatorAperture: public RectEllipseAperture
{
public:
	/**
	 *  CollimatorAperture default constructor
	 *  @param[in] w collimator FULL entrance width
	 *  @param[in] h collimator FULL entrance height
	 *  @param[in] t collimator tilt angle
	 *  @param[in] length collimator length
	 *  @param[in] x_offset_entry horizontal offset
	 *  @param[in] y_offset_entry vertical offset
	 */
	CollimatorAperture(double w, double h, double t, double length, double x_offset_entry = 0.0, double y_offset_entry =
		0.0);

	/**
	 *  function to set collimator aperture FULL entrance width
	 *  @param[in] collimator aperture entrance width
	 */
	void SetEntranceWidth(double);

	/**
	 *  function to set collimator aperture FULL entrance height
	 *  @param[in] collimator aperture entrance height
	 */
	void SetEntranceHeight(double);

	/**
	 *  function to set collimator aperture FULL exit width
	 *  @param[in] collimator aperture exit width
	 */
	void SetExitWidth(double);

	/**
	 *  function to set collimator aperture FULL exit height
	 *  @param[in] collimator aperture exit height
	 */
	void SetExitHeight(double);

	/**
	 *  function to set collimator aperture exit horizontal offset
	 *  @param[in] collimator aperture exit horizontal offset
	 */
	void SetExitXOffset(double);

	/**
	 *  function to set collimator aperture exit vertical offset
	 *  @param[in] collimator aperture exit vertical offset
	 */
	void SetExitYOffset(double);

	/**
	 *  function to get collimator aperture FULL entrance width
	 *  @return collimator aperture entrance width
	 */
	double GetFullEntranceWidth() const;

	/**
	 *  function to get collimator aperture FULL entrance height
	 *  @return collimator aperture entrance height
	 */
	double GetFullEntranceHeight() const;

	/**
	 *  function to get collimator aperture FULL exit width
	 *  @return collimator aperture exit width
	 */
	double GetFullExitWidth() const;

	/**
	 *  function to get collimator aperture FULL exit height
	 *  @return collimator aperture exit height
	 */
	double GetFullExitHeight() const;

	/**
	 *  function to get collimator aperture horizontal offset at entrance
	 *  @return horizontal offset at entrance
	 */
	double GetEntranceXOffset() const;

	/**
	 *  function to get collimator aperture vertical offset at entrance
	 *  @return vertical offset at entrance
	 */
	double GetEntranceYOffset() const;

	/**
	 *  function to get collimator aperture horizontal offset at exit
	 *  @param[in] horizontal offset at exit
	 */
	double GetExitXOffset() const;

	/**
	 *  function to get collimator aperture vertical offset at exit
	 *  @return vertical offset at exit
	 */
	double GetExitYOffset() const;

	/**
	 *  function to get collimator tilt
	 *  @return collimator tilt angle
	 */
	double GetCollimatorTilt() const;

	/**
	 *  function to get collimator length
	 *  @return collimator length
	 */
	double GetCollimatorLength() const;

	/**
	 *  CollimatorAperture virtual function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	virtual bool CheckWithinApertureBoundaries(double x, double y, double z) const;

protected:
	double alpha;
	double CollimatorLength;
	double x_offset_entry, y_offset_entry;
	double x_offset_exit, y_offset_exit;
	double w_entrance, h_entrance;
	double w_exit, h_exit;
	double cosalpha;
	double sinalpha;

};

class UnalignedCollimatorAperture: public CollimatorAperture
{
public:

	/**
	 * Default constructor
	 */
	UnalignedCollimatorAperture(double w, double h, double t, double length, double x_offset_entry = 0.0, double
		y_offset_entry = 0.0);

	/**
	 *  UnalignedCollimatorAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;
};

class CollimatorApertureWithErrors: public CollimatorAperture
{
	double ApertureError;

	/**
	 *  CollimatorApertureWithErrors override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;
};

class UnalignedCollimatorApertureWithErrors: public UnalignedCollimatorAperture
{
	double ApertureError;

	/**
	 *  UnalignedCollimatorApertureWithErrors override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;
};

class OneSidedUnalignedCollimatorAperture: public CollimatorAperture
{
public:
	/**
	 * Default constructor
	 */
	OneSidedUnalignedCollimatorAperture(double w, double h, double t, double length, double x_offset_entry = 0.0, double
		y_offset_entry = 0.0, bool side = 1);

	/**
	 *  UnalignedCollimatorApertureWithErrors override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;
	bool JawSide;

	/**
	 * function to set Jaw side
	 * @param[in] bool jaw side - 0 outer, 1 inner - check??
	 */
	void SetJawSide(bool);

	/**
	 * function to get Jaw side
	 * @return bool collimator jaw side - 0 outer, 1 inner - check??
	 */
	bool GetJawSide();
};

#endif
