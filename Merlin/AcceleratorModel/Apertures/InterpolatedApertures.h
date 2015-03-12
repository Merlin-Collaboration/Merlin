#ifndef _InterpolatedApertures_h_
#define _InterpolatedApertures_h_ 1

#include "SimpleApertures.h"
#include <vector>
#include <string>

class InterpolatedAperture
{
public:

	struct ap
	{
		double s;
		double ap1;
		double ap2;
		double ap3;
		double ap4;
	};
	InterpolatedAperture(){};
	InterpolatedAperture(std::vector<ap> ApertureListInput) : ApertureList(ApertureListInput) {};

	std::vector<ap> GetApertureList() const
	{
		return ApertureList;
	}
	void AddS(double s)
	{
		ApertureEntry.s = s;
	}
	void AddAp1(double ap1)
	{
		ApertureEntry.ap1 = ap1;
	}
	void AddAp2(double ap2)
	{
		ApertureEntry.ap2 = ap2;
	}
	void AddAp3(double ap3)
	{
		ApertureEntry.ap3 = ap3;
	}
	void AddAp4(double ap4)
	{
		ApertureEntry.ap4 = ap4;
	}
	void AddEntry()
	{
		ApertureList.push_back(ApertureEntry);
	}
//protected:

	ap ApertureEntry;
	std::vector<ap> ApertureList;
};

//RectEllipse Aperture
class InterpolatedRectEllipseAperture : public InterpolatedAperture, public Aperture
{
public:
	//Constructor
	InterpolatedRectEllipseAperture (std::vector<InterpolatedAperture::ap> ElementApertureList):InterpolatedAperture(ElementApertureList)
	{
		//MerlinProfile::AddProcess("APERTURE");
		Print = false;
		ApertureMaterial = NULL;
	}

	//Returns true if the point (x,y,z) is within the aperture.
	virtual bool PointInside (double x, double y, double z) const;
    
	//TODO
	//Returns true if the point p is within the aperture.
	//bool PointInside (const Point3D& p) const;

	//Returns the radius to the aperture at location z and angle phi.
	virtual double GetRadiusAt (double phi, double z) const;
	//double GetFullHeight () const;
	//double GetFullWidth () const;

	virtual std::string GetApertureType() const;

	void EnablePrint();
	bool Print;
	virtual void printout(std::ostream& out) const;

/*
	void SetRectHalfWidth(double);
	void SetRectHalfHeight(double);
	void SetEllipseHalfHorizontal(double);
	void SetEllipseHalfVertical(double);

//private:
	double rect_half_width;
	double rect_half_height;
	double ellipse_half_horizontal;
	double ellipse_half_vertical;
*/
};

inline string InterpolatedRectEllipseAperture::GetApertureType() const
{
	return "INTERPOLATEDRECTELLIPSE";
}

class InterpolatedCircularAperture : public InterpolatedAperture, public Aperture
{
public:
	InterpolatedCircularAperture (std::vector<ap> ElementApertureList):InterpolatedAperture(ElementApertureList) {};

	double GetRadius () const;
	double GetDiameter () const;
	void SetRadius (double r);
	void SetDiameter (double d);

	//Returns true if the point (x,y,z) is within the aperture.
	virtual bool PointInside (double x, double y, double z) const;

	//Returns the radius.
	virtual double GetRadiusAt (double phi, double z) const;

	virtual std::string GetApertureType() const;
	virtual void printout(std::ostream& out) const;
private:
	double r2;
};

inline string InterpolatedCircularAperture::GetApertureType() const
{
	return "INTERPOLATEDCIRCULAR";
}
#endif
