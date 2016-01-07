#include "AcceleratorModel/Apertures/CollimatorAperture.h"
#include "MADInterface/MADInterface.h"
#include "Random/RandomNG.h"


/**********************************************************************
*
*	A collimtor jaw, aligned to the beam orbit and beta function changes
* 	This does NOT have jaw flatness errors
*
**********************************************************************/
CollimatorAperture::CollimatorAperture(double w,double h, double t, Material* m, double length, double x_off, double y_off)
	:RectangularAperture(w,h), alpha(t), CollimatorLength(length), x_offset_entry(x_off), y_offset_entry(y_off),cosalpha(cos(-t)),sinalpha(sin(-t))
{
	SetMaterial(m);
	x_offset_exit = 0;
	y_offset_exit = 0;
	w_exit = 0;
	h_exit = 0;
	//cout << cosalpha << "\t" << sinalpha << endl;
}

//Checks if particle is in or outside a defined aperture
inline bool CollimatorAperture::PointInside(double x,double y,double z) const
{
	/*
		if(errors)
		{
			//Should use/adjust GetFullWidth
			x += aperture_error * ((pow(z,2)/jaw_length)-z);
			y += aperture_error * ((pow(z,2)/jaw_length)-z);
		}
	*/

	/*
	We need to calculate several variables:
	1: The x,y offsets at the z position of the particle.
	2: The width and hight of the jaw at the z position of the particle.

	Lets start with the x and y offsets;
	*/

	//y = mx + c
	//m = dy/dx

	double x_off = (z * ( x_offset_entry - x_offset_exit ) / CollimatorLength) - x_offset_entry;
	double y_off = (z * ( y_offset_entry - y_offset_exit ) / CollimatorLength) - y_offset_entry;


	//These will give the jaw width and heights to be used. * 0.5 to convert to half width.
	double x_jaw = (z * ( w_exit - GetFullWidth() )  / CollimatorLength) + GetFullWidth();
	double y_jaw = (z * ( h_exit - GetFullHeight() ) / CollimatorLength) + GetFullHeight();

	double x1 = ((x+x_off) * cosalpha) - ((y+y_off) * sinalpha);
	double y1 = ((x+x_off) * sinalpha) + ((y+y_off) * cosalpha);

	//output everything
	//~ if(! (fabs(x1) < (x_jaw/2) && fabs(y1) < (y_jaw/2)) ){
	//~ cout << "\nCollAp: z = " << z << "\t\t CollimatorLength = " << CollimatorLength << endl;
	//~ cout << "x_off = " << x_off << "\t\t x_offset_entry = " << x_offset_entry << "\t\t x_offset_entry " << x_offset_entry << endl;
	//~ //cout << "y_off = " << y_off << "\t\t y_offset_entry = " << y_offset_entry << "\t\t y_offset_entry " << y_offset_entry << endl;
	//~ cout << "x_jaw = " << x_jaw << z << "\t\t w_exit = " << w_exit << "\t\t Full Width = " << GetFullWidth() << endl;
	//~ //cout << "y_jaw = " << y_jaw << z << "\t\t h_exit = " << h_exit << "\t\t GetFullHeight() " << GetFullHeight() << endl;
	//~ cout << "x = " << x << "\t\t x_off = " << x_off << "\t\t cosalpha = " << cosalpha << endl;
	//~ cout << "x1 = " << x1 << "\t\t x_jaw/2 = " << x_jaw/2 << endl;
	//~ }

	return fabs(x1) < (x_jaw/2) && fabs(y1) < (y_jaw/2);
	/*	ofstream outf;
		outf.precision(18);
		outf.open("/samdata2/serlucam/Output/jaw.dat",fstream::app);
		outf <<  z << "\t" <<  x_off << "\t"  <<  y_off << "\t" <<  x_jaw << "\t" << y_jaw << "\t" << alpha << "\t" << x1 << endl;
		outf.close();
	*/
}


//Sets the jaw width at the exit of the collimator
void CollimatorAperture::SetExitWidth(double width)
{
	w_exit = width;
}

//Sets the jaw height at the exit of the collimator
void CollimatorAperture::SetExitHeight(double height)
{
	h_exit = height;
}

//Sets the x (horizontal) orbit offset at the exit of the collimator
void CollimatorAperture::SetExitXOffset(double x)
{
	x_offset_exit = x;
}

//Sets the y (vertical) orbit offset at the exit of the collimator
void CollimatorAperture::SetExitYOffset(double y)
{
	y_offset_exit = y;
}
/*
//Sets the collimator length
void CollimatorAperture::SetCollimatorLength(double length)
{
	CollimatorLength = length;
}

void CollimatorAperture::SetJawLength(double length)
{
	jaw_length = length;
}
*/

double CollimatorAperture::GetFullEntranceHeight() const
{
	return GetFullHeight();
}

double CollimatorAperture::GetFullEntranceWidth() const
{
	return GetFullWidth();
}

double CollimatorAperture::GetFullExitHeight() const
{
	return h_exit;
}

double CollimatorAperture::GetFullExitWidth() const
{
	return w_exit;
}

double CollimatorAperture::GetEntranceXOffset() const
{
	return x_offset_entry;
}

double CollimatorAperture::GetEntranceYOffset() const
{
	return y_offset_entry;
}

double CollimatorAperture::GetExitXOffset() const
{
	return x_offset_exit;
}

double CollimatorAperture::GetExitYOffset() const
{
	return y_offset_exit;
}

double CollimatorAperture::GetCollimatorTilt() const
{
	return alpha;
}


/**********************************************************************
*
*	A collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This does NOT have jaw flatness errors
*	Still aligned to the beam size! Flat jaws, parallel to the beam pipe,
*	but touching the beam envelope on each side
*	This is the LHC configuration
*
**********************************************************************/
UnalignedCollimatorAperture::UnalignedCollimatorAperture(double w,double h, double t, Material* m, double length, double x_off, double y_off)
	:CollimatorAperture(w,h,t,m,length,x_off,y_off)
//:RectangularAperture(w,h), alpha(t), CollimatorLength(length), x_offset(x_off), y_offset(y_off)
{
	SetMaterial(m);
}

inline bool UnalignedCollimatorAperture::PointInside(double x,double y,double z) const
{
	double x1 = ((x-x_offset_entry) * cosalpha) - ((y-y_offset_entry) * sinalpha);
	double y1 = ((x-x_offset_entry) * sinalpha) + ((y-y_offset_entry) * cosalpha);
//	cout << "aperture checking" << endl;

	if(fabs(x1) < GetFullWidth()/2)
	{
		//	cout << "passed x: " << x1 << "\t" << GetFullWidth()/2 << endl;
	}
	else
	{
		//cout << "failed x: " << x1 << "\t" << x << "\t" <<x_offset_entry << "\t" <<  GetFullWidth()/2 << "\t" << alpha << endl;
	}
	if(fabs(y1) < GetFullHeight()/2)
	{
		//	cout << "passed y: " << y1 << "\t" << GetFullHeight()/2 << endl;
	}
	else
	{
		//cout << "failed y: " << y1 << "\t" << GetFullHeight()/2 << endl;
	}

	return fabs(x1) < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
}


/**********************************************************************
*
*	A collimtor jaw, aligned to the beam orbit and beta function changes
* 	This has jaw flatness errors
*
**********************************************************************/

inline bool CollimatorApertureWithErrors::PointInside(double x,double y,double z) const
{
	double x_off = (z * ( x_offset_entry - x_offset_exit ) / CollimatorLength) - x_offset_entry;
	double y_off = (z * ( y_offset_entry - y_offset_exit ) / CollimatorLength) - y_offset_entry;

	//These will give the jaw width and heights to be used. * 0.5 to convert to half width.
	double x_jaw = (z * ( w_exit - GetFullWidth() )  / CollimatorLength) + GetFullWidth();
	double y_jaw = (z * ( h_exit - GetFullHeight() ) / CollimatorLength) + GetFullHeight();

	double x1 = ((x+x_off) * cosalpha) - ((y+y_off) * sinalpha);
	double y1 = ((x+x_off) * sinalpha) + ((y+y_off) * cosalpha);
	return fabs(x1) < (x_jaw/2) && fabs(y1) < (y_jaw/2);
	x1 += ApertureError * ((pow(z,2)/jaw_length)-z);
	y1 += ApertureError * ((pow(z,2)/jaw_length)-z);

	return fabs(x1) < (x_jaw/2) && fabs(y1) < (y_jaw/2);
}

/**********************************************************************
*
*	A collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This has jaw flatness errors
*
**********************************************************************/
/*
inline bool UnalignedCollimatorApertureWithErrors::PointInside(double x,double y,double z) const
{
	double x1 = ((x+x_offset_entry) * cos(-alpha)) - ((y+y_offset_entry) * sin(-alpha));
	double y1 = ((x+x_offset_entry) * sin(-alpha)) + ((y+y_offset_entry) * cos(-alpha));
	x1 += ApertureError * ((pow(z,2)/jaw_length)-z);
	y1 += ApertureError * ((pow(z,2)/jaw_length)-z);

	return fabs(x1) < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
}
*/

/**********************************************************************
*
*	A one sided collimtor jaw, unaligned to the beam orbit or beta function changes
* 	This does NOT have jaw flatness errors
*	Still aligned to the beam size! Flat jaws, parallel to the beam pipe,
*	but touching the beam envelope on each side
*	This is the LHC configuration
*
**********************************************************************/
OneSidedUnalignedCollimatorAperture::OneSidedUnalignedCollimatorAperture(double w,double h, double t, Material* m, double length, double x_off, double y_off)
	:CollimatorAperture(w,h,t,m,length,x_off,y_off),PositiveSide(true)
//:RectangularAperture(w,h), alpha(t), CollimatorLength(length), x_offset(x_off), y_offset(y_off)
{
	SetMaterial(m);
}

inline bool OneSidedUnalignedCollimatorAperture::PointInside(double x,double y,double z) const
{
	double x1 = ((x-x_offset_entry) * cosalpha) - ((y-y_offset_entry) * sinalpha);
	double y1 = ((x-x_offset_entry) * sinalpha) + ((y-y_offset_entry) * cosalpha);
	//cout << "aperture checking" << endl;

	if(PositiveSide)
	{
		return x1 < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
	}
	else
	{
		return (-x1) < GetFullWidth()/2 && fabs(y1) < GetFullHeight()/2;
	}
}

void OneSidedUnalignedCollimatorAperture::SetJawSide(bool side)
{
	PositiveSide = side;
}
