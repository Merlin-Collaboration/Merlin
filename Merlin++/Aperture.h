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
#include <map>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "DataTable.h"

class Aperture
{
public:

	/**
	 *  Aperture default constructor
	 *  @param[in] type the aperture typename string
	 *  @param[in] s the longitudinal location coordinate
	 *  @param[in] aper1 the horizontal rectangular aperture parameter
	 *  @param[in] aper2 the vertical rectangular aperture parameter
	 *  @param[in] aper3 the horizontal elliptical aperture parameter
	 *  @param[in] aper4 the vertical elliptical aperture parameter
	 */
	Aperture();
	Aperture(std::string type, double s, double aper1);
	Aperture(std::string type, double s, double aper1, double aper2);
	Aperture(std::string type, double s, double aper1, double aper2, double aper3, double aper4);

	/**
	 *  Virtual destructor
	 */
	virtual ~Aperture();

	/**
	 *	Pure virtual function/interface for CheckWithinApertureBoundaries
	 *  @param[in] x the x location coordinate of the particle
	 *  @param[in] y the y location coordinate of the particle
	 *  @param[in] z the z location coordinate of the particle
	 *  @return bool flag confirming particle is within boundaries
	 */
	virtual bool CheckWithinApertureBoundaries(double x, double y, double z) const = 0;

	/**
	 *	Pure virtual function/interface for getting aperture type
	 *	@return string of aperture typename
	 */
	virtual std::string GetType() = 0;

	/**
	 *	Function/interface for setting aperture type
	 *	@return string of aperture typename
	 */
	void SetType(std::string type)
	{
		apType = type;
	}

	/**
	 *  function to set longitudinal location coordinate
	 *  @param[in] s he longitudinal location coordinate
	 */
	void SetLatticeLocation(double s)
	{
		latticelocation = s;
	}

	/**
	 *  function to set longitudinal location coordinate
	 *  @param[in] s he longitudinal location coordinate
	 */
	double GetLatticeLocation() const
	{
		return latticelocation;
	}

	/**
	 *  function to allow printing of all read apertures
	 *  @param[in] out The stream to print
	 */
	virtual void printout(std::ostream& out) const;

protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	std::string apType;
	double latticelocation;
};

class CircularAperture: public Aperture
{
public:
	/**
	 *  CircularAperture simple constructor
	 */
	CircularAperture(double radius);

	/**
	 *  CircularAperture override base constructor
	 */
	CircularAperture(std::string type, double s, double radius);

	/**
	 *  CircularAperture override of Aperture member function getType()
	 */
	std::string GetType();

	/**
	 * function to get the radius aperture parameter
	 *  @return the radius aperture parameter
	 */
	double GetRadius()
	{
		return radius;
	}

	/**
	 *  CircularAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new CircularAperture instance - only called by ApertureFactory::GetInstance class
	 *  @return constructed Aperture pointer of assigned type CircularAperture
	 */
	//static Aperture* GetInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
	static Aperture* GetInstance(DataTableRow);

	virtual void printout(std::ostream& out) const;
protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	double radius;
	double radius_sq;
};

class RectangularAperture: public Aperture
{
public:

	/**
	 *  RectangularAperture simple constructor
	 */
	RectangularAperture(double rectHalfX, double rectHalfY);

	/**
	 *  RectangularAperture override default constructor
	 */
	RectangularAperture(std::string type, double s, double rectHalfX, double rectHalfY);

	/**
	 *  RectangularAperture override of Aperture member function getType()
	 */
	std::string GetType();

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfX()
	{
		return rectHalfX;
	}

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfY()
	{
		return rectHalfY;
	}

	virtual void printout(std::ostream& out) const;

	/**
	 *  RectangularAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new RectangularAperture instance - only called by ApertureFactory::GetInstance class
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	//static Aperture* GetInstance(string type, double s, double aper1, double aper2, double aper3, double aper4);
	static Aperture* GetInstance(DataTableRow);

protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	double rectHalfX;
	double rectHalfY;
	double minDim;
	double maxDim;
};

class EllipticalAperture: public Aperture
{
public:

	/**
	 *  EllipticalAperture override base constructor
	 */
	EllipticalAperture(std::string type, double s, double ellipHalfX, double ellipHalfY);

	/**
	 *  EllipticalAperture override of Aperture member function getType()
	 */
	std::string GetType();

	/**
	 *  function to get the horizontal elliptical aperture parameter
	 *  @return the horizontal elliptical aperture parameter
	 */
	double GetEllipHalfX()
	{
		return ellipHalfX;
	}

	/**
	 *  function to get the vertical elliptical aperture parameter
	 *  @return the vertical elliptical aperture parameter
	 */
	double GetEllipHalfY()
	{
		return ellipHalfY;
	}
	/**
	 *  EllipticalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new EllipticalAperture instance - only called by ApertureFactory::GetInstance class
	 *  @return constructed Aperture pointer of assigned type EllipticalAperture
	 */
	static Aperture* GetInstance(DataTableRow);

	virtual void printout(std::ostream& out) const;
protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	double ellipHalfX;
	double ellipHalfY;
	double minDim;
	double maxDim;
	double ellipHalfX_sq;
	double ellipHalfY_sq;
	double ellipHalfX_sq_div_ellipHalfY_sq;
};

class RectEllipseAperture: public Aperture
{
public:

	/**
	 *  RectEllipseAperture default constructor for collimator aperture inheritence
	 */
	RectEllipseAperture();

	/**
	 *  RectEllipseAperture override default constructor
	 */
	RectEllipseAperture(std::string type, double s, double rectHalfX, double rectHalfY, double ellipHalfX, double
		ellipHalfY);

	/**
	 *  RectEllipseAperture override of Aperture member function getType()
	 */
	std::string GetType();

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfX()
	{
		return rectHalfX;
	}

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfY()
	{
		return rectHalfY;
	}

	/**
	 *  function to get the horizontal elliptical aperture parameter
	 *  @return the horizontal elliptical aperture parameter
	 */
	double GetEllipHalfX()
	{
		return ellipHalfX;
	}

	/**
	 *  function to get the vertical elliptical aperture parameter
	 *  @return the vertical elliptical aperture parameter
	 */
	double GetEllipHalfY()
	{
		return ellipHalfY;
	}

	/**
	 *  RectEllipseAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new RectEllipseAperture instance - only called by ApertureFactory::GetInstance class
	 *  @return constructed Aperture pointer of assigned type RectEllipseAperture
	 */
	static Aperture* GetInstance(DataTableRow);

	struct ap
	{
		std::string ApType;
		double s;
		double ap1;
		double ap2;
		double ap3;
		double ap4;

	};

	ap ApertureEntry;

	std::vector<ap> ApertureList;

	virtual void printout(std::ostream& out) const;
protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	double rectHalfX;
	double rectHalfY;
	double ellipHalfX;
	double ellipHalfY;
	double minDim;
	double maxDim;
	double minRectDim;
	double maxRectDim;
	double minEllipDim;
	double maxEllipDim;
	double ellipHalfX_sq;
	double ellipHalfY_sq;
	double ellipHalfX_sq_div_ellipHalfY_sq;
};

class OctagonalAperture: public Aperture
{
public:

	/**
	 *  OctagonalAperture override default constructor
	 *  Note: For octagonal apertures, aper3 and aper4 are angles, not dimensions
	 *  see (slide 6): https://indico.cern.ch/event/379692/contributions/1804923/subcontributions/156446/attachments/757501/1039118/2105-03-18_HSS_meeting_rev.pdf
	 *
	 */
	OctagonalAperture(std::string type, double s, double rectHalfX, double rectHalfY, double angle1, double angle2);

	/**
	 *  OctagonalAperture override of Aperture member function getType()
	 */
	std::string GetType();

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfX()
	{
		return rectHalfX;
	}

	/**
	 * function to get the vertical rectangular aperture parameter
	 *  @return the vertical rectangular aperture parameter
	 */
	double GetRectHalfY()
	{
		return rectHalfY;
	}

	/**
	 *  function to get the angle1 octagonal aperture parameter
	 *  @return the angle1 octagonal aperture parameter
	 */
	double GetAngle1()
	{
		return angle1;
	}

	/**
	 *  function to get the angle2 octagonal aperture parameter
	 *  @return the angle2 octagonal aperture parameter
	 */
	double GetAngle2()
	{
		return angle2;
	}

	/**
	 *  OctagonalAperture override of Aperture member function CheckWithinApertureBoundaries()
	 *  @param[in] x x-coord of particle
	 *  @param[in] y y-coord of particle
	 *  @param[in] z z-coord of particle
	 *  @return true/false flag
	 */
	bool CheckWithinApertureBoundaries(double x, double y, double z) const;

	/**
	 *  get new OctagonalAperture instance - only called by ApertureFactory::GetInstance class
	 *  @return constructed Aperture pointer of assigned type OctagonalAperture
	 */
	static Aperture* GetInstance(DataTableRow);

	virtual void printout(std::ostream& out) const;
protected:

	/**
	 *  protected aperture parameters accessible only via class get/set functions
	 */
	double rectHalfX;
	double rectHalfY;
	double minDim;
	double maxDim;
	double angle1;
	double angle2;
	double const1;
	double const2;
	double const3;
};

/**
 * typedef of default constructor to getAperture member function pointer
 */
typedef Aperture* (*getAperture)(DataTableRow);

class ApertureFactory
{
public:
	/**
	 * define map of typename string to type-specific constructor member functions
	 */
	static std::map<std::string, getAperture> ApertureTypes;

	/**
	 * gets instance of input typename-specific Aperture, checks type against ApertureFactoryInitializer list
	 * @return constructed Aperture pointer of assigned type
	 */
	static Aperture* GetInstance(DataTableRow);
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
