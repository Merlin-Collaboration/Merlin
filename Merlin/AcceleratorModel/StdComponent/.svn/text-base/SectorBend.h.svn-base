/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (2000)
 * 
 * file Merlin\AcceleratorModel\StdComponent\SectorBend.h
 * last modified 04/04/01 15:24:09
 */

/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 *
 * Copyright (c) 2000 by The Merlin Collaboration.  
 * ALL RIGHTS RESERVED. 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */


#ifndef SectorBend_h
#define SectorBend_h 1

#include "merlin_config.h"


// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"
// ArcGeometry
#include "AcceleratorModel/StdGeometry/ArcGeometry.h"
// ArcMultipoleField
#include "AcceleratorModel/StdComponent/ArcMultipoleField.h"

class ComponentTracker;




//	Representation of a SectorBend magnet. By default, the
//	magnet is hard edged, with the pole-faces perpendicular
//	to the geometry. Pole faces can be rotated about the
//	local y-axis, and can also have a constant radius
//	curvature whioch gives rise to a sextupole component
//	(see TRANSPORT or MAD for details). As well as a
//	vertical B field, sector bends can have an additional
//	quadrupole and sextupole field component specified.


class SectorBend : public ArcMultipoleField
{
public:
    //	POD containing information for the pole face of a sector
    //	bend.

    class PoleFace
    {
    public:
        PoleFace (double angle = 0, double f_int = 0, double hg = 0);

        // Data Members for Class Attributes

        //	Pole face rotation.
        double rot;

        //	fringe field integral.
        double fint;

        //	half gap magnet.
        double hgap;

    protected:
    private:
    private:
    };

    //	Pair struct for containing pointers to the entrance and
    //	exit pole face attributes.


    class PoleFaceInfo
    {
    public:
        PoleFaceInfo ()
                : entrance(0),exit(0)
        {
        }

        PoleFaceInfo (const PoleFaceInfo& rhs)
        {
            Copy(rhs);
        }

        ~PoleFaceInfo ()
        {
            Clear();
        }


        PoleFaceInfo& operator = (const PoleFaceInfo& rhs)
        {
            Clear();
            Copy(rhs);
            return *this;
        }

        void SetInfo (PoleFace* e1, PoleFace* e2)
        {
            Clear();
            entrance=e1;
            exit=e2;
        }

        void SetInfo (PoleFace* e1)
        {
            Clear();
            entrance=exit=e1;
        }

        // Data Members for Class Attributes

        PoleFace* entrance;

        PoleFace* exit;

    protected:
    private:
    private:

        void Clear ()
        {
            if(entrance) delete entrance;
            if(exit && exit!=entrance) delete exit;
        }

        void Copy (const PoleFaceInfo& rhs);

    };

public:
    //	Constructor taking the id of the dipole, the length
    //	(meter) and curvature (1/meter) of the dipole geometry,
    //	and the vertical magnetic field in Tesla.
    SectorBend (const string& id, double len, double h, double b0);


    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    virtual const string& GetType () const;

    //	Rotate the bend by 180 degrees about the vertical
    //	access. Causes the sign of the curvature to change, as
    //	well as the field to be rotated (via a call to Multipole
    //	Magnet::rotate()).
    virtual void RotateY180 ();

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    virtual ModelElement* Copy () const;

    const PoleFaceInfo& GetPoleFaceInfo () const
    {
        return pfInfo;
    }

    void SetPoleFaceInfo (PoleFace* entr, PoleFace* exit);

    void SetPoleFaceInfo (PoleFace* pf);

    //	Returns the matched momentum  for the current field and
    //	geometry for a particle of charge q/e. Note that the
    //	returned momentum can be negative.
    double GetMatchedMomentum (double q) const;

    //	Returns the main dipole field in Tesla.
    double GetB0 () const
    {
        return GetField().GetFieldScale();
    }

    //	Returns the quadrupole field component  in Tesla/meter.
    double GetB1 () const;

    //	Sets the main dipole field (in Tesla).
    void SetB0 (double By)
    {
        GetField().SetFieldScale(By);
    }

    //	Sets the quadrupole field component for the dipole
    //	(mixed function magnet). B1 is defined as the quadrupole
    //	gradient in Tesla/meter.
    void SetB1 (double b1);

    // Data Members for Class Attributes

    //	Unique index for an Accelerator component.
    static const int ID;

	// The followind field access function added for
	// compatability with other magnets
	void SetFieldStrength(double b) { SetB0(b); }
	double GetFieldStrength() const { return GetB0(); }

protected:
private:
    // Data Members for Associations

    PoleFaceInfo pfInfo;

private:
};



#endif
