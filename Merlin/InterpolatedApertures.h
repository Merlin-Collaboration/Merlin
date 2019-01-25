/*
 * InterpolatedAperture.h
 *
 *  Created on: 31 Oct 2017
 *      Author: scott
 */

#ifndef INTERPOLATEDAPERTURE_H_
#define INTERPOLATEDAPERTURE_H_

#include "Aperture.h"

using namespace std;

/**
 * Interpolated apertures utilize all four parameters and are inherently of rectellipse geometry
 */

class InterpolatedRectEllipseAperture: public Aperture
{
public:

	/**
	 *  InterpolatedRectEllipseAperture default constructor
	 *  @param[in] apVec vector of relevant aperture pointers to be interpolated
	 */
	InterpolatedRectEllipseAperture(vector<Aperture*> apVec);

	/**
	 *  Virtual destructor
	 */
	virtual ~InterpolatedRectEllipseAperture();

	/**
	 *  function to get aperture typename string
	 *  @return aperture typename string
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
	 *  get new InterpolatedRectEllipseAperture instance - only called by ApertureFactory::getInstance class
	 *  @param[in] vector of Aperture
	 *  @return constructed Aperture pointer of assigned type RectangularAperture
	 */
	static Aperture* getInstance(vector<Aperture*>);

	vector<Aperture*> ElementApertures;
};

typedef Aperture* (*getInterpolator)(vector<Aperture*>);

class InterpolatorFactory
{
public:
	static map<string, getInterpolator> interpolatorTypes;
	Aperture* getInstance(vector<Aperture*>);
};

class InterpolatorFactoryInitializer
{
	static InterpolatorFactoryInitializer init;
public:
	InterpolatorFactoryInitializer();
};

#endif /* INTERPOLATEDAPERTURE_H_ */
