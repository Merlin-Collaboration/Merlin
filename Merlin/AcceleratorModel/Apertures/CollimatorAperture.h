#ifndef _CollimatorAperture_h_
#define _CollimatorAperture_h_

#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "Collimators/MaterialDatabase.h"
#include "Collimators/Material.h"
#include <iostream>

/**********************************************************************
*
*	A collimtor jaw, aligned to the beam orbit and beta function changes
* 	This does NOT have jaw flatness errors
*
**********************************************************************/

class CollimatorAperture: public RectangularAperture
{
protected:
	double alpha;
	double CollimatorLength;
	double jaw_length;
	double x_offset_entry,y_offset_entry;

//Add jaw parameters at exit as well
	double x_offset_exit,y_offset_exit;
	double w_exit,h_exit;
	double cosalpha;
	double sinalpha;
//bool errors;
//double e1,e2,e3,e4,e5,e6;
//double aperture_error;
//void EnableErrors(bool);
//void SetErrors(double, double);


public:
//~ double alpha;
	CollimatorAperture(double w,double h, double t, Material* m, double length, double x_offset_entry=0.0, double y_offset_entry=0.0);

//void SetJawLength(double);

	void SetExitWidth(double);	//Horizontal
	void SetExitHeight(double);	//Vertical
	void SetExitXOffset(double);	//Horizontal
	void SetExitYOffset(double);	//Vertical

	double GetFullEntranceHeight() const;
	double GetFullEntranceWidth() const;

	double GetFullExitHeight() const;
	double GetFullExitWidth() const;

	double GetEntranceXOffset() const;
	double GetEntranceYOffset() const;

	double GetExitXOffset() const;
	double GetExitYOffset() const;

	double GetCollimatorTilt() const;

//Also need to know the collimator length for interpolation
//void SetCollimatorLength(double);

	virtual bool PointInside(double x,double y,double z) const;
};


/**********************************************************************
*
*	A collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This does NOT have jaw flatness errors
*
**********************************************************************/

class UnalignedCollimatorAperture: public CollimatorAperture
{
	/*
	protected:
	double alpha;
	double CollimatorLength;
	double jaw_length;
	double x_offset,y_offset;
	*/
public:
	UnalignedCollimatorAperture(double w,double h, double t, Material* m, double length, double x_offset_entry=0.0, double y_offset_entry=0.0);

	bool PointInside(double x,double y,double z) const;
};

/**********************************************************************
*
*	A collimtor jaw, aligned to the beam orbit and beta function changes
* 	This has jaw flatness errors
*
**********************************************************************/

class CollimatorApertureWithErrors: public CollimatorAperture
{
	double ApertureError;
	bool PointInside(double x,double y,double z) const;
};

/**********************************************************************
*
*	A collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This has jaw flatness errors
*
**********************************************************************/

class UnalignedCollimatorApertureWithErrors: public UnalignedCollimatorAperture
{
	double ApertureError;
	bool PointInside(double x,double y,double z) const;
};


/**********************************************************************
*
*	A collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This does NOT have jaw flatness errors
*
**********************************************************************/

class OneSidedUnalignedCollimatorAperture: public CollimatorAperture
{
	/*
	protected:
	double alpha;
	double CollimatorLength;
	double jaw_length;
	double x_offset,y_offset;
	*/
public:
	OneSidedUnalignedCollimatorAperture(double w,double h, double t, Material* m, double length, double x_offset_entry=0.0, double y_offset_entry=0.0);

	bool PointInside(double x,double y,double z) const;
	bool PositiveSide;
	void SetJawSide(bool);
};

#endif
