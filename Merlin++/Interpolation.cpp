/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Interpolation.h"
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;

// Method used for arbitrarily space data
//
class ArbSpacedData: public Interpolation::Method
{
public:

	struct Data
	{
		double x, y;
		Data(double x1 = 0, double y1 = 0) :
			x(x1), y(y1)
		{
		}

	};

	int lastseenat;

	ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals, int order);
	ArbSpacedData(const double* xvals, const double* yvals, int n, int order);
	virtual double ValueAt(double, double* err = 0x0);

private:
	vector<Data> itsData;
	vector<vector<double> > Q;
};

// Method used for equally space data
//
class EqualSpacedData: public Interpolation::Method
{
public:

	EqualSpacedData(const vector<double> yv, double xm, double delta, int order) :
		yvals(yv), xmin(xm), xmax(xm + (yv.size() - 1) * delta), dx(delta)
	{
		itsOrder = order;
		assert(dx > 0);
	}

	double ValueAt(double x, double*err = 0x0);

private:
	vector<double> yvals;
	double xmin;
	double xmax;
	double dx;
};

// Function  used by sort() for ArbSpacedData::Data
//
inline bool sort_x(const ArbSpacedData::Data& d1, const ArbSpacedData::Data& d2)
{
	return d1.x < d2.x;
}

// Class ArbSpacedData implementation
//

ArbSpacedData::ArbSpacedData(const vector<double>& xvals, const vector<double>& yvals, int order) :
	itsData()
{
	assert(xvals.size() == yvals.size());
	itsData.reserve(xvals.size());
	for(size_t i = 0; i < xvals.size(); i++)
	{
		itsData.push_back(Data(xvals[i], yvals[i]));
	}
	sort(itsData.begin(), itsData.end(), sort_x);

	itsOrder = order;
	int n = xvals.size();
	if(itsOrder > n - 1)
	{
		itsOrder = n - 1;
		cout << " Interpolation with only " << n << " points so requested order reduced from " << order << " to " <<
			itsOrder << endl;
	}

	//  Neville's method
	//  First, for  n points define n constants
	//  then define n-1 linear forms between adjacent pairs
	//  then use these to give n-2 quadratics that go through triplets
	//  etcetera

	vector<vector<double> > P;
	P.reserve(n);

	for(int i = 0; i < n; i++)
	{
		vector<double> v;
		v.reserve(n);
		v.push_back(itsData[i].y);
		P.push_back(v);
	}

	// this is the iteration, from old P to new Q
	for(int k = 1; k <= itsOrder; k++)
	{
		Q.clear();
		Q.reserve(n - k);
		for(int i = 0; i < n - k; i++)
		{
			Q.push_back({0}); // needs this to get size() right later
			Q[i].clear();
			Q[i].reserve(k + 1);
			for(int j = 0; j < k; j++)
			{
				Q[i].push_back((-itsData[i + k].x * P[i][j] + itsData[i].x * P[i + 1][j]) / (itsData[i].x - itsData[i
					+ k].x));
			}
			Q[i].push_back(0);
			for(int j = 1; j <= k; j++)
			{
				Q[i][j] += (P[i][j - 1] - P[i + 1][j - 1]) / (itsData[i].x - itsData[i + k].x);
			}
		}

		P.clear();         // replace old P with copy of new Q for use in next iteration
		P.reserve(n - k);
		for(int i = 0; i < n - k; i++)
		{
			P.push_back(Q[i]);
		}
	}
	lastseenat = 0;  // used for finding segment: assumes repeated nearby calls
	/*
	   if(itsOrder<2) return;
	   // DEBUG HERE
	   cout<<n<<" points\n";
	   cout<<" order "<<itsOrder<<endl;
	   cout<<" Qsize "<<Q.size()<<endl;
	   cout<<" Psize "<<P.size()<<endl;
	   std::ofstream f("jim.txt");
	    double last=xvals[xvals.size()-1];

	    double stepp=(last-xvals[0])/100;
	    for(double xx=xvals[0];xx<last;xx+=stepp){
	            f<<xx;
	            for(uint ii=0;ii<Q.size();ii++) f<<" "<<evaluatepoly(xx,Q[ii]);
	            f<<endl;
	            }

	   f.close();
	 */
}
ArbSpacedData::ArbSpacedData(const double* xvals, const double* yvals, int n, int order) :
	itsData()
{ // version for old style arrays
	vector<double> X(xvals, xvals + n);
	vector<double> Y(yvals, yvals + n);
	ArbSpacedData(X, Y, order);
}

double ArbSpacedData::ValueAt(double x, double* err)
{
	if(x < itsData.front().x || x > itsData.back().x)
	{
		throw Interpolation::BadRange(x, FloatRange(itsData.front().x, itsData.back().x));
	}

	int n = itsData.size();

	if(itsOrder > 1)
	{
		// this covers interpolation more complicated than linear
		// search is inefficient first time but fast for many close calls
		// n.b.  lastseenat is a data member and therefore preserved between calls
		while(x > itsData[lastseenat].x)
			lastseenat += 1;
		while(x < itsData[lastseenat].x)
			lastseenat -= 1;
		int lo = max(0, lastseenat - itsOrder);
		int Qsize = Q.size();
		int hi = min(lastseenat, Qsize - 1);
		int iwant = max(0, lastseenat - (itsOrder + 1) / 2 - 1);
		iwant = min(iwant, Qsize - 2);
		// cout<<" lastseen lo hi iwant "<<lastseenat<<" "<<lo<<" "<<hi<<" "<<iwant<<endl;
		//for(int i=lo;i<=hi;i++) cout<<i<<" "<<evaluatepoly(x,Q[i])<<endl;
		double v1 = evaluatepoly(x, Q[iwant]);
		double v2 = evaluatepoly(x, Q[iwant + 1]);
		//cout<<" is "<<v1 <<" "<<v2<<endl;
		if(err)
			*err = abs(v1 - v2);  // give an error estimate iff one is asked for
		return (v1 + v2) / 2;
	}
	// use linear interpolation to return value
	// locate segment by binary search
	size_t ju = itsData.size(), jl = 0, jm;

	while((ju - jl) > 1)
	{
		jm = (ju + jl) >> 1;
		if(x > itsData[jm].x)
		{
			jl = jm;
		}
		else
		{
			ju = jm;
		}
	}
	double m = (itsData[jl + 1].y - itsData[jl].y) / (itsData[jl + 1].x - itsData[jl].x);
	return itsData[jl].y + (x - itsData[jl].x) * m;
}

// Class EqualSpacedData implementation
//

double EqualSpacedData::ValueAt(double x, double* err)
{
	// note that we use extrapolation here if x is out of range
	size_t n;
	if(x < xmin)
	{
		n = 0;
	}
	else if(x > xmax)
	{
		n = yvals.size() - 2;
	}
	else
	{
		n = (x - xmin) / dx;
	}

	double m = (yvals[n + 1] - yvals[n]) / dx;
	double x0 = xmin + dx * n;
	return yvals[n] + m * (x - x0);
}
// exception
Interpolation::BadRange::BadRange(double x, const FloatRange& r) :
	MerlinException("BadRange")
{
	ostringstream buf;
	buf << x << " not in interpolation range (" << r.lower << "," << r.upper << ')';
	msg = buf.str();
}

Interpolation::Interpolation(const vector<double>& yvals, double xmin, double dx, int order) :
	itsMethod(new EqualSpacedData(yvals, xmin, dx, order))
{
}

Interpolation::Interpolation(const double* xv, const double* yv, int n, int order) :
	itsMethod(new ArbSpacedData(xv, yv, n, order))
{
}
Interpolation::Interpolation(const std::vector<double>& xv, const vector<double>& yv, int order) :
	itsMethod(new ArbSpacedData(xv, yv, order))
{
}

Interpolation::~Interpolation()
{
	if(itsMethod)
	{
		delete itsMethod;
	}
}
