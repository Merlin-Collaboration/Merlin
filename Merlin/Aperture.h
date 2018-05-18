/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Aperture_h
#define Aperture_h 1

#include "merlin_config.h"
#include "Space3D.h"
#include "Material.h"
#include <iostream>
#include <string>

/**
* Represents the cross section of the vacuum pipe or other
* collimating aperture.
*/
class Aperture
{
public:

	/**
	* Constructor
	* @param[in] m The material to attach to this aperture.
	*/
	Aperture(Material* m = nullptr) : ApertureMaterial(m) {}

	/**
	* Destructor
	*/
	virtual ~Aperture ();

	/**
	* Returns true if the point (x,y,z) is within the aperture.
	*
	* @param[in] x The x coordinate of the particle
	* @param[in] y The y coordinate of the particle
	* @param[in] z The z coordinate of the particle
	* @retval true If the specified point is within the Aperture
	* @retval false.
	*/
	virtual bool PointInside (double x, double y, double z) const = 0;

	/**
	* Returns true if the point p is within the aperture.
	*
	* @param[in] p 3D point reference to check
	* @retval true If the specified point is within the Aperture
	* @retval false
	*/
	bool PointInside (const Point3D& p) const;

	/**
	* Returns the radius to the aperture at location z and angle phi.
	* @param[in] phi The angle.
	* @param[in] z The z position of the particle.
	* @return A double containing the aperture radius.
	*/
	virtual double GetRadiusAt (double phi, double z) const = 0;

	/**
	* Returns the type of the aperture.
	* @return A string containing the type of the aperture.
	*/
	virtual std::string GetApertureType() const = 0;

	/**
	* Gets the material associated with this aperture.
	* @return A pointer to a Material class attached to this aperture.
	*/
	Material* GetMaterial() const;

	/**
	* Sets the material associated with this aperture.
	* @param[in] m The material to attach to this aperture.
	*/
	void SetMaterial(Material* m);

	/**
	* Prints out the Aperture parameters to a specified stream.
	* @param[out] out The output stream to use.
	*/
	virtual void printout(std::ostream& out) const;

protected:

	/**
	* The Material of this aperture.
	*/
	Material* ApertureMaterial;
};

std::ostream& operator<< (std::ostream& out, const Aperture& ap);

/**
* See the MAD users guide for how these apertures are defined.
* (current as of V5.02.07)
* http://madx.web.cern.ch/madx/releases/last-dev/madxuguide.pdf
* "Physical Aperture: Aperture definition"
*
* Interpolated in this case is where one type joins another - future internal usage, not a MAD-X type.
*/
typedef enum
{
	NONE,
	UNKNOWN,
	CIRCLE,			//Supported
	RECTANGLE,		//Supported
	ELLIPSE,		//Supported
	RECTCIRCLE,
	LHCSCREEN,		//Supported as RECTELLIPSE
	RECTELLIPSE,	//Supported
	RACETRACK,
	OCTAGON,
	INTERPOLATED
} ApertureClass;

#endif
