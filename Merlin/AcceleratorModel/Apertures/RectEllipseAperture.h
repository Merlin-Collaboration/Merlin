#ifndef _RectEllipseAperture_h_
#define _RectEllipseAperture_h_

#include "merlin_config.h"
#include <cmath>
#include "AcceleratorModel/Aperture.h"

/*
* Rectellipse Aperture
* This class defines the LHC beam screen
*/
class RectEllipseAperture : public Aperture
{
public:
        //Constructor
        RectEllipseAperture (double rhw, double rhh, double ehh, double ehv)
        : Aperture(),RectHalfWidth(rhw),RectHalfHeight(rhh),
		EllipseHalfHorizontal(ehh),EllipseHalfVertical(ehv),
		EHH2(ehh*ehh),HV((ehh*ehh)/(ehv*ehv)){};
        //: Aperture(),RectHalfWidth(rhw),RectHalfHeight(rhh),EllipseHalfHorizontal(ehh),EllipseHalfVertical(ehv),EHH2(ehh*ehh),HV((ehh*ehh)/(ehv*ehv)) {};

        //Returns true if the point (x,y,z) is within the
        //aperture. The z coordinate is ignored.
        bool PointInside (double x, double y, double z) const;

        //Returns the radius to the aperture at location z and angle phi.
        double GetRadiusAt (double phi, double z) const;
        std::string GetApertureType() const;
		virtual void printout(std::ostream& out) const;

protected:
	const double RectHalfWidth;
	const double RectHalfHeight;
	const double EllipseHalfHorizontal;
	const double EllipseHalfVertical;

	/*
	* The squares of the above numbers
	*/
	/*
	const double RHW2;
	const double RHH2;
	const double EHH2;
	const double EHV2;
	*/

	/*
	* The inverse of the above numbers
	* IRHW2 = 1 / (RHW*RHW), etc
	*/
	/*
	const double IRHW2;
	const double IRHH2;
	const double IEHH2;
	const double IEHV2;
	*/

	//Just using the below two save a multiply operation
	//pow(EllipseHalfHorizontal,2)
	const double EHH2;
	//pow(EllipseHalfHorizontal,2) / pow(EllipseHalfVertical,2)
	const double HV;
};

#endif
