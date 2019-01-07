/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef APERTURE_H_
#define APERTURE_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

using namespace std;

class ApertureAbstract
{
public:

	/**
	 *  Abstract constructor
	 */
	ApertureAbstract();

	/**
	 *  Virtual destructor
	 */
	virtual ~ApertureAbstract();

	/**
	 *	Pure virtual function/interface for CheckWithinApertureBoundaries
	 *  @param[in] x the x location coordinate of the particle
	 *  @param[in] y the y location coordinate of the particle
	 *  @param[in] z the z location coordinate of the particle
	 *  @return bool flag confirming particle is within boundaries
	 */
	virtual bool CheckWithinApertureBoundaries(double x, double y, double z) const = 0;

	/**
	 *	Pure virtual function/interface for getting particle type
	 *	@return string of aperture typename
	 */
	virtual string getType() = 0;
};

class Aperture: public ApertureAbstract
{
public:

	/**
	 *  constructor
	 */
	Aperture();

	/**
	 *  Aperture default constructor
	 *  @param[in] type the aperture typename string
	 *  @param[in] s the longitudinal location coordinate
	 *  @param[in] aper1 the horizontal rectangular aperture parameter
	 *  @param[in] aper2 the vertical rectangular aperture parameter
	 *  @param[in] aper3 the horizontal elliptical aperture parameter
	 *  @param[in] aper4 the vertical elliptical aperture parameter
	 */
	Aperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  function to set aperture type variable
	 *  @param[in] apType the typename string you wish to set
	 */
	void setApertureType(string aperType)
	{
		apType = aperType;
	}

	/**
	 *  function to set longitudinal location coordinate
	 *  @param[in] s he longitudinal location coordinate
	 */
	void setSlongitudinal(double s)
	{
		s_longitudinal = s;
	}

	/**
	 *  function to calculate all common aperture parameters
	 *  @param[in] aper1 the horizontal rectangular aperture parameter
	 *  @param[in] aper2 the vertical rectangular aperture parameter
	 *  @param[in] aper3 the horizontal elliptical aperture parameter
	 *  @param[in] aper4 the vertical elliptical aperture parameter
	 */
	void CalcApertureParams(double aper1, double aper2, double aper3, double aper4)
	{
		minRectDim = fmin(aper1, aper2);
		minEllipDim = fmin(aper3, aper4);
		minDim = min({aper1, aper2, aper3, aper4});
		maxRectDim = fmax(aper1, aper2);
		maxEllipDim = fmax(aper3, aper4);
		maxDim = max({aper1, aper2, aper3, aper4});
		ellipHalfHeight2 = aper4 * aper4;
		ellipHalfWidth2 = aper3 * aper3;
		ellipHalfWidth2overEllipHalfHeight2 = ellipHalfWidth2 / ellipHalfHeight2;
	}

	/**
	 *  function to set the horizontal rectangular aperture parameter
	 *  @param[in] aper1 the horizontal rectangular aperture parameter
	 */
	void setRectHalfWidth(double aper1)
	{
		rectHalfWidth = aper1;
		minRectDim = fmin(rectHalfWidth, rectHalfHeight);
		maxRectDim = fmax(rectHalfWidth, rectHalfHeight);
		minDim = min({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		maxDim = max({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
	}

	/**
	 *  function to set the vertical rectangular aperture parameter
	 *  @param[in] aper2 the vertical rectangular aperture parameter
	 */
	void setRectHalfHeight(double aper2)
	{
		rectHalfHeight = aper2;
		minRectDim = fmin(rectHalfWidth, rectHalfHeight);
		maxRectDim = fmax(rectHalfWidth, rectHalfHeight);
		minDim = min({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		maxDim = max({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
	}

	/**
	 *  function to set the horizontal elliptical aperture parameter
	 *  @param[in] aper3 the horizontal elliptical aperture parameter
	 */
	void setEllipHalfWidth(double aper3)
	{
		ellipHalfWidth = aper3;
		minEllipDim = fmin(ellipHalfWidth, ellipHalfHeight);
		maxEllipDim = fmax(ellipHalfWidth, ellipHalfHeight);
		minDim = min({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		maxDim = max({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		ellipHalfWidth2 = ellipHalfWidth * ellipHalfWidth;
		ellipHalfWidth2overEllipHalfHeight2 = ellipHalfWidth2 / ellipHalfHeight2;
	}

	/**
	 *  function to set the vertical elliptical aperture parameter
	 *  @param[in] aper4 the vertical elliptical aperture parameter
	 */
	void setEllipHalfHeight(double aper4)
	{
		ellipHalfHeight = aper4;
		minEllipDim = fmin(ellipHalfWidth, ellipHalfHeight);
		maxEllipDim = fmax(ellipHalfWidth, ellipHalfHeight);
		minDim = min({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		maxDim = max({rectHalfWidth, rectHalfHeight, ellipHalfWidth, ellipHalfHeight});
		ellipHalfWidth2 = ellipHalfHeight * ellipHalfHeight;
		ellipHalfWidth2overEllipHalfHeight2 = ellipHalfWidth2 / ellipHalfHeight2;
	}

	/**
	 * function to set min dimension of aperture for - used in collimators
	 * @param[in] minDim the minimum aperture dimension
	 */

	void setMinDim(double mindim)
	{
		minDim = mindim;
	}

	/**
	 *  function to get aperture typename string
	 *  @return the aperture typename string
	 */
	string getApertureType()
	{
		return apType;
	}

	/**
	 *  function to get aperture typename string from factory instance
	 *  @return the aperture typename string from factory instance
	 */
	string getType()
	{
		return apType;
	}

	/**
	 *  function to get the longitudinal location coordinate
	 *  @return the longitudinal location coordinate
	 */
	double getSlongitudinal()
	{
		return s_longitudinal;
	}

	/**
	 *  function to get the horizontal rectangular aperture parameter
	 *  @return the horizontal rectangular aperture parameter
	 */
	double getRectHalfWidth()
	{
		return rectHalfWidth;
	}

	/**
	 *  function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double getRectHalfHeight()
	{
		return rectHalfHeight;
	}

	/**
	 *  function to get the horizontal elliptical aperture parameter
	 *  @return the horizontal elliptical aperture parameter
	 */
	double getEllipHalfWidth()
	{
		return ellipHalfWidth;
	}

	/**
	 *  function to get the vertical elliptical aperture parameter
	 *  @return the vertical elliptical aperture parameter
	 */
	double getEllipHalfHeight()
	{
		return ellipHalfHeight;
	}

	/**
	 *  function to allow printing of all read apertures
	 *  @param[in] out The stream to print
	 */
	void printout(std::ostream& out) const
	{
		out << apType;
	}

protected:
	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */

	string apType;
	double s_longitudinal;
	double rectHalfWidth;
	double rectHalfHeight;
	double ellipHalfWidth;
	double ellipHalfHeight;
	double minRectDim;
	double minEllipDim;
	double minDim;
	double maxRectDim;
	double maxEllipDim;
	double maxDim;
	double ellipHalfHeight2;
	double ellipHalfWidth2;
	double ellipHalfWidth2overEllipHalfHeight2;
};

class CircularAperture: public Aperture
{
public:
	/**
	 *  CircularAperture override base constructor
	 */
	CircularAperture(double aper3);

	/**
	 *  CircularAperture override default constructor
	 */
	CircularAperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  CircularAperture override of Aperture member function getType()
	 */
	string getType();

	/**
	 *  CircularAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new CircularAperture instance - only called by ApertureFactory::getInstance class
	 *  @return constructed Aperture pointer of assigned type CircularAperture
	 */
	static Aperture* getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
};

class RectangularAperture: public Aperture
{
public:
	/**
	 *  RectangularAperture override base constructor
	 */
	RectangularAperture(double aper1, double aper2);

	/**
	 *  RectangularAperture override default constructor
	 */
	RectangularAperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  RectangularAperture override of Aperture member function getType()
	 */
	string getType();

	/**
	 *  RectangularAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new RectangularAperture instance - only called by ApertureFactory::getInstance class
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
};

class EllipticalAperture: public Aperture
{

public:

	/**
	 *  EllipticalAperture override base constructor
	 */
	EllipticalAperture(double aper3, double aper4);

	/**
	 *  EllipticalAperture override default constructor
	 */
	EllipticalAperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  EllipticalAperture override of Aperture member function getType()
	 */
	string getType();

	/**
	 *  EllipticalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new EllipticalAperture instance - only called by ApertureFactory::getInstance class
	 *  @return constructed Aperture pointer of assigned type EllipticalAperture
	 */
	static Aperture* getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
};

class RectEllipseAperture: public Aperture
{

public:

	/**
	 *  RectEllipseAperture empty constructor for collimator override
	 */
	RectEllipseAperture();
	/**
	 *  RectEllipseAperture override base constructor
	 */
	RectEllipseAperture(double aper1, double aper2, double aper3, double aper4);

	/**
	 *  RectEllipseAperture override default constructor
	 */
	RectEllipseAperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  RectEllipseAperture override of Aperture member function getType()
	 */
	string getType();

	/**
	 *  RectEllipseAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new RectEllipseAperture instance - only called by ApertureFactory::getInstance class
	 *  @return constructed Aperture pointer of assigned type RectEllipseAperture
	 */
	static Aperture* getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
};

class OctagonalAperture: public Aperture
{

public:

	/**
	 *  OctagonalAperture override base constructor
	 */
	OctagonalAperture();

	/**
	 *  OctagonalAperture override default constructor
	 *  Note: For octagonal apertures, aper3 and aper4 are angles, not dimensions
	 *  see (slide 6): https://indico.cern.ch/event/379692/contributions/1804923/subcontributions/156446/attachments/757501/1039118/2105-03-18_HSS_meeting_rev.pdf
	 *
	 */
	OctagonalAperture(string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  OctagonalAperture override of Aperture member function getType()
	 */
	string getType();

	/**
	 *  OctagonalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new OctagonalAperture instance - only called by ApertureFactory::getInstance class
	 *  @return constructed Aperture pointer of assigned type OctagonalAperture
	 */
	static Aperture* getInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);

protected:
	/**
	 * Octagon angles, see MAD site: for details
	 */

	double angle1;
	double angle2;

	/**
	 * Aperture calculation constants
	 */

	double const1;
	double const2;
	double const3;
};

/**
 * typedef of default constructor to getAperture member function pointer
 */
typedef Aperture* (*getAperture)(string type, double s, double aper1, double aper2, double aper3, double aper4);

class ApertureFactory
{
public:
	/**
	 * define map of typename string to type-specific constructor member functions
	 */
	static map<string, getAperture> apertureTypes;

	/**
	 * gets instance of input typename-specific Aperture, checks type against ApertureFactoryInitializer list
	 * @return constructed Aperture pointer of assigned type
	 */
	Aperture* getInstance(string type = "RECTELLIPSE", double s = 0.0, double aper1 = 0.0, double aper2 = 0.0, double
		aper3 = 0.0, double aper4 = 0.0);
};

class ApertureFactoryInitializer
{
	static ApertureFactoryInitializer init;
public:

	/**
	 * ApertureFactoryInitializer constructor
	 * contains list of viable Aperture types for use by the ApertureFactory
	 */
	ApertureFactoryInitializer();
};

#endif /* APERTURE_H_ */
