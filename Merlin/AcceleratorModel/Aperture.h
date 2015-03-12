/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef Aperture_h
#define Aperture_h 1

#include "merlin_config.h"
#include "EuclideanGeometry/Space3D.h"
#include "Collimators/Material.h"
#include <iostream>
#include <string>

/**
* Represents the cross section of the vacuum pipe or other
* collimating aperture.
*/
class Aperture
{
public:
	Aperture(Material* m = NULL) : ApertureMaterial(m){}

	virtual ~Aperture ();

	/**
	* Returns true if the point (x,y,z) is within the aperture.
	* @param x The x coordinate of the particle
	* @param y The y coordinate of the particle
	* @param z The z coordinate of the particle
	*/
	virtual bool PointInside (double x, double y, double z) const = 0;

	/**
	* Returns true if the point p is within the aperture.
	* @param p 3D point reference to check
	*/
	bool PointInside (const Point3D& p) const;

	/**
	* Returns the radius to the aperture at location z and angle phi.
	* @param phi The angle
	* @param z The z posision of the particle
	*/
	virtual double GetRadiusAt (double phi, double z) const = 0;

	/**
	* Returns the type of the aperture.
	*/
	virtual string GetApertureType() const = 0;

	Material* GetMaterial() const {return ApertureMaterial;}
	void SetMaterial(Material* m){ApertureMaterial = m ;}
	virtual void printout(std::ostream& out) const;

protected:
	Material* ApertureMaterial;
};

inline Aperture::~Aperture () {}

inline bool Aperture::PointInside (const Point3D& p) const
{
	return PointInside(p.x,p.y,p.z);
}


std::ostream& operator<< (std::ostream& out, const Aperture& ap);

#endif
