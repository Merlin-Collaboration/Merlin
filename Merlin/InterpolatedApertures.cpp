/*
 * InterpolatedAperture.cpp
 *
 *  Created on: 31 Oct 2017
 *      Author: scott
 */

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "InterpolatedApertures.h"

using namespace std;

string InterpolatedRectEllipseAperture::getType()
{
	return "RECTELLIPSEinterpolated";
}

InterpolatedRectEllipseAperture::InterpolatedRectEllipseAperture(vector<Aperture*> apVec) :
	ElementApertures(apVec)
{

}

InterpolatedRectEllipseAperture::~InterpolatedRectEllipseAperture()
{

}

bool InterpolatedRectEllipseAperture::CheckWithinApertureBoundaries(double x, double y, double z) const
{
	Aperture* apBack;
	Aperture* apFront;
	double ax = fabs(x);
	double ay = fabs(y);

	for(size_t n = 1; n < ElementApertures.size(); n++)
	{
		if(ElementApertures[n]->getSlongitudinal() >= z)
		{
			apFront = ElementApertures[n];
			apBack = ElementApertures[n - 1];
			break;
		}
	}
	if(ax + ay
		< min(
			{ apFront->getRectHalfWidth(), apFront->getRectHalfHeight(), apFront->getEllipHalfWidth(),
			  apFront->getEllipHalfHeight(), apBack->getRectHalfWidth(), apBack->getRectHalfHeight(),
			  apBack->getEllipHalfWidth(), apBack->getEllipHalfHeight() }))
		return true;

	double delta_s = apFront->getSlongitudinal() - apBack->getSlongitudinal();
	double rectHalfWidthAtZ = apFront->getRectHalfWidth()
		- (((apFront->getRectHalfWidth() - apBack->getRectHalfWidth()) / delta_s)
		* (apFront->getSlongitudinal() - z));
	double rectHalfHeithAtZ = apFront->getRectHalfHeight()
		- (((apFront->getRectHalfHeight() - apBack->getRectHalfHeight()) / delta_s)
		* (apFront->getSlongitudinal() - z));
	double ellipHalfWidthAtZ = apFront->getEllipHalfWidth()
		- (((apFront->getEllipHalfWidth() - apBack->getEllipHalfWidth()) / delta_s)
		* (apFront->getSlongitudinal() - z));
	double ellipHalfHeightAtZ = apFront->getEllipHalfHeight()
		- (((apFront->getEllipHalfHeight() - apBack->getEllipHalfHeight()) / delta_s)
		* (apFront->getSlongitudinal() - z));

	if(ax + ay < min(
			{ rectHalfWidthAtZ, rectHalfHeithAtZ, ellipHalfWidthAtZ, ellipHalfHeightAtZ }))
		return true;
	if(((x * x) / (ellipHalfWidthAtZ * ellipHalfWidthAtZ)) + ((y * y) / (ellipHalfHeightAtZ * ellipHalfHeightAtZ)) > 1)
		return false;
	if(ax > rectHalfWidthAtZ || ay > rectHalfHeithAtZ)
		return false;
	else
		return true;
}

Aperture* InterpolatedRectEllipseAperture::getInstance(vector<Aperture*> apVec)
{
	return new InterpolatedRectEllipseAperture(apVec);
}

Aperture* InterpolatorFactory::getInstance(vector<Aperture*> apVec)
{
	map<string, getInterpolator>::iterator itr = interpolatorTypes.find(apVec[0]->getType());
	if(itr != interpolatorTypes.end())
	{
		return (*itr->second)(apVec);
	}
	return NULL;
}

InterpolatorFactoryInitializer::InterpolatorFactoryInitializer()
{
//	InterpolatorFactory::apertureTypes["CIRCLE"] = &newCircularAperture::getInstance;
//	InterpolatorFactory::apertureTypes["RECTANGLE"] = &newRectangularAperture::getInstance;
//	InterpolatorFactory::apertureTypes["ELLIPSE"] = &newEllipticalAperture::getInstance;
	InterpolatorFactory::interpolatorTypes["RECTELLIPSE"] = &InterpolatedRectEllipseAperture::getInstance;
//	InterpolatorFactory::apertureTypes["OCTAGON"] = &newOctagonalAperture::getInstance;
}

map<string, getInterpolator> InterpolatorFactory::interpolatorTypes;
InterpolatorFactoryInitializer InterpolatorFactoryInitializer::init;
