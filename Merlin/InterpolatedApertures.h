/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef INTERPOLATEDAPERTURE_H_
#define INTERPOLATEDAPERTURE_H_

#include "Aperture.h"
#include "DataTable.h"

using namespace std;

/**
 * Interpolated apertures utilize all four parameters and are inherently of rectellipse geometry
 */

class InterpolatedAperture: public Aperture
{
public:

	/**
	 *  InterpolatedAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedAperture(DataTable ApertureDataTable);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
	 */
	string GetType();

	/**
	 *  InterpolatedAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new InterpolatedAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* GetInstance(DataTable);

	DataTable AperturesToInterpolate;

	struct apStruct
	{
		string ApType;
		double s;
		double ap1;
		double ap2;
		double ap3;
		double ap4;

	};

	apStruct ApEntry;

	vector<apStruct> ApList;

	void ConvertToStruct(DataTable);
};

class InterpolatedRectEllipseAperture: public InterpolatedAperture
{
public:
	/**
	 *  InterpolatedRectEllipseAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedRectEllipseAperture(DataTable ApertureDataTable);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedRectEllipseAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
	 */
	string GetType();

	/**
	 *  InterpolatedRectEllipseAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new InterpolatedRectEllipseAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* GetInstance(DataTable);
};

class InterpolatedCircularAperture: public InterpolatedAperture
{
public:
	/**
	 *  InterpolatedCircularAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedCircularAperture(DataTable ApertureDataTable);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedCircularAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
	 */
	string GetType();

	/**
	 *  InterpolatedCircularAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new InterpolatedCircularAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* GetInstance(DataTable);
};

class InterpolatedEllipticalAperture: public InterpolatedAperture
{
public:
	/**
	 *  InterpolatedEllipticalAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedEllipticalAperture(DataTable ApertureDataTable);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedEllipticalAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
	 */
	string GetType();

	/**
	 *  InterpolatedEllipticalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new InterpolatedEllipticalAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* GetInstance(DataTable);
};

class InterpolatedOctagonalAperture: public InterpolatedAperture
{
public:
	/**
	 *  InterpolatedOctagonalAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedOctagonalAperture(DataTable ApertureDataTable);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedOctagonalAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
	 */
	string GetType();

	/**
	 *  InterpolatedOctagonalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new InterpolatedOctagonalAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* GetInstance(DataTable);
};

typedef Aperture* (*getInterpolator)(DataTable);

class InterpolatorFactory
{
public:
	static map<string, getInterpolator> interpolatorTypes;
	Aperture* GetInstance(DataTable);
};

class InterpolatorFactoryInitializer
{
	static InterpolatorFactoryInitializer init;
public:
	InterpolatorFactoryInitializer();
};

#endif /* INTERPOLATEDAPERTURE_H_ */
