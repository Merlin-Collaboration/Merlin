/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SectorBend_h
#define SectorBend_h 1

#include "merlin_config.h"

#include "TemplateComponents.h"
#include "MultipoleField.h"
#include "ArcGeometry.h"
#include "ArcMultipoleField.h"

class ComponentTracker;

/**
 *	Representation of a SectorBend magnet. By default, the
 *	magnet is hard edged, with the pole-faces perpendicular
 *	to the geometry. Pole faces can be rotated about the
 *	local y-axis, and can also have a constant radius
 *	curvature which gives rise to a sextupole component
 *	(see TRANSPORT or MAD for details). As well as a
 *	vertical B field, sector bends can have an additional
 *	quadrupole and sextupole field component specified.
 */
class SectorBend: public ArcMultipoleField
{
public:

	/**
	 *	POD containing information for the pole face of a sector
	 *	bend.
	 */
	class PoleFace
	{
	public:
		PoleFace(double angle = 0, double f_int = 0, double hg = 0, double type = 1);

		/**
		 *	Pole face rotation -- data member for class attributes
		 */
		double rot;

		/**
		 *	fringe field integral -- data member for class attributes
		 */
		double fint;

		/**
		 *	half gap magnet -- data member for class attributes
		 */
		double hgap;

		/**
		 *	Pole face type -- data member for class attributes
		 */
		double type;

	protected:
	private:
	private:
	};

	/**
	 *	Pair struct for containing pointers to the entrance and
	 *	exit pole face attributes.
	 */
	class PoleFaceInfo
	{
	public:
		PoleFaceInfo() :
			entrance(nullptr), exit(nullptr)
		{
		}

		PoleFaceInfo(const PoleFaceInfo& rhs)
		{
			Copy(rhs);
		}

		~PoleFaceInfo()
		{
			Clear();
		}

		PoleFaceInfo& operator =(const PoleFaceInfo& rhs)
		{
			Clear();
			Copy(rhs);
			return *this;
		}

		void SetInfo(PoleFace* e1, PoleFace* e2)
		{
			Clear();
			entrance = e1;
			exit = e2;
			entrance->type = 1;
			exit->type = 0;
		}

		void SetInfo(PoleFace* e1)
		{
			Clear();
			entrance = exit = e1;
			entrance->type = 1;
			exit->type = 0;
		}

		/**
		 * Data member for class attributes
		 */
		PoleFace* entrance;

		/**
		 * Data member for class attributes
		 */
		PoleFace* exit;

	protected:
	private:
	private:

		void Clear()
		{
			if(entrance)
			{
				delete entrance;
			}
			if(exit && exit != entrance)
			{
				delete exit;
			}
		}

		void Copy(const PoleFaceInfo& rhs);

	};

public:
	/**
	 *	Constructor taking the id of the dipole, the length
	 *	(meter) and curvature (1/meter) of the dipole geometry,
	 *	and the vertical magnetic field in Tesla.
	 */
	SectorBend(const string& id, double len, double h, double b0);

	/**
	 *	Returns the unique index for this class of accelerator
	 *	components.
	 */
	virtual int GetIndex() const;

	virtual const string& GetType() const;

	/**
	 *	Rotate the bend by 180 degrees about the vertical
	 *	access. Causes the sign of the curvature to change, as
	 *	well as the field to be rotated (via a call to Multipole
	 *	Magnet::rotate()).
	 */
	virtual void RotateY180();

	/**
	 *	Primary tracking interface. Prepares the specified
	 *	Tracker object for tracking this component.
	 */
	virtual void PrepareTracker(ComponentTracker& aTracker);

	virtual ModelElement* Copy() const;

	const PoleFaceInfo& GetPoleFaceInfo() const
	{
		return pfInfo;
	}

	void SetPoleFaceInfo(PoleFace* entr, PoleFace* exit);

	void SetPoleFaceInfo(PoleFace* pf);

	/**
	 *	Returns the matched momentum  for the current field and
	 *	geometry for a particle of charge q/e. Note that the
	 *	returned momentum can be negative.
	 *
	 *	@return Matched momentum for current field and geometry
	 */
	double GetMatchedMomentum(double q) const;

	/**
	 *	Returns the main dipole field in Tesla.
	 *	@return Main dipole field
	 */
	double GetB0() const
	{
		return GetField().GetFieldScale();
	}

	/**
	 *	Returns the quadrupole field component  in Tesla/meter.
	 *	@return Quadrupole field component
	 */
	double GetB1() const;

	/**
	 *	Sets the main dipole field (in Tesla).
	 */
	void SetB0(double By)
	{
		GetField().SetFieldScale(By);
	}

	/**
	 *	Sets the quadrupole field component for the dipole
	 *	(mixed function magnet). B1 is defined as the quadrupole
	 *	gradient in Tesla/meter.
	 */
	void SetB1(double b1);

	// Data Members for Class Attributes

	/**
	 *	Unique index for an Accelerator component.
	 */
	static const int ID;

	/**
	 * The following field access function added for
	 * compatibility with other magnets
	 */
	void SetFieldStrength(double b)
	{
		SetB0(b);
	}
	double GetFieldStrength() const
	{
		return GetB0();
	}

protected:
private:
	// Data Members for Associations

	PoleFaceInfo pfInfo;

private:
};

#endif
