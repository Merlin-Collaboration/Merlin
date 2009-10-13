/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/07 09:14:12 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef Interpolation_h
#define Interpolation_h 1

#include "merlin_config.h"
#include "Exception/MerlinException.h"
#include "NumericalUtils/Range.h"
#include <vector>

// class Interpolation
// An interpolation functor which interpolates values from a data table.
// Currenly only linear interpolation is assumed.

class Interpolation {
public:

    // exception
class BadRange : public MerlinException {
    public:
        BadRange(double x, const FloatRange& r);
        double value;
        FloatRange valid_range;
    };

    // Implementation method for interpolation
    class Method {
    public:
        virtual ~Method() {}
        virtual double ValueAt(double x) const throw (BadRange) =0;
    };

    // Interpolation of equally spaced data points
    Interpolation(const std::vector<double>& yvals, double xmin, double dx);

    // Interpolation of arbitrary spaced data points
    Interpolation(const std::vector<double>& xvals, const std::vector<double>& yvals);

    ~Interpolation();

    double operator()(double x) const throw (BadRange) {
        return itsMethod->ValueAt(x);
    }

private:

    Method* itsMethod;
};

#endif

