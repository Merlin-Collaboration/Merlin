/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef EquilibriumDistribution_h
#define EquilibriumDistribution_h 1

// A.Wolski 19/11/01

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/StdComponent/SectorBend.h"

class IntegrateEigenvector
{
public:
    double Integral(ComplexVector& Ek, SectorBend* sb, double p0);

protected:
    ComplexVector ek;
    Complex k;
    double h, tanE1;

    void polint(double xa[], double ya[], int n, double x, double& y, double& dy);
    double trapzd(double a, double b, int n);
    double qromb(double a, double b);
    virtual double func(double s) = 0;
};

class IntegrateZeroGradient : public IntegrateEigenvector
{
public:
    IntegrateZeroGradient();
protected:
    virtual double func(double s);
};

class IntegrateWithGradient : public IntegrateEigenvector
{
public:
    IntegrateWithGradient();
protected:
    virtual double func(double s);
};

class EquilibriumDistribution
{
public:
    EquilibriumDistribution(AcceleratorModel* aModel, double refMomentum);
    double Emittance(int n);
    double DampingTime(int n);
    double DampingConstant(int n);
    double Tune(int n);
    double SynchronousTime();
    double BeamMoment(int i, int j, int ncpt=0);

    void CalculateDampingConstants();
    void CalculateEmittance();

private:
    AcceleratorModel* theModel;
    double p0;
    double dampingConstant[3];
    double tune[3];
    double emittance[3];
    double synchronousTime;

};

#endif

