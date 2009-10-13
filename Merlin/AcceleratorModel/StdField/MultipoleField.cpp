/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iterator>
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"

using namespace std;

namespace {

double factorial(int n)
{
    double f=n--;
    while(n>1) f*=n--;
    return (f==0) ? 1 : f;  //Conditional added by A.Wolski 1/10/01
}

}; // end anonymous namespace


MultipoleField::MultipoleField (int np, double bn, double an, double r0)
        :B0(),expansion(np+1,Complex(0,0))
{
    Complex z(bn,an);
    B0=abs(z);
    expansion[np]=z/B0;
    B0/=pow(r0,np);
}

MultipoleField::MultipoleField (int np, double bn, double r0, bool skew)
        : B0(),expansion(np+1,Complex(0,0))
{
    expansion[np] = skew ? Complex(0,1) : Complex(1,0);
    B0=bn/pow(r0,np);
}

MultipoleField::MultipoleField (int np, double bn, bool skew)
        :B0(bn/factorial(np)),expansion(np+1,Complex(0,0))
{
    expansion[np] = skew ? Complex(0,1) : Complex(1,0);
}

double MultipoleField::GetFieldScale () const
{
    return B0;
}

void MultipoleField::SetFieldScale (double scale)
{
    B0=scale;
}

Complex MultipoleField::GetKn (int np, double rigidity) const
{
    assert(rigidity!=0);
    return static_cast<size_t>(np)<expansion.size() ? (B0/rigidity)*factorial(np)*expansion[np] : Complex(0);
}

Complex MultipoleField::GetField2D (double x, double y, int exclude) const
{
    if(IsNullField() || (++exclude)>static_cast<int>(expansion.size()))
        return Complex(0);

    const Complex z0(x,y);

    Complex B(0,0);
    Complex z(1);

    for(int n=0; n<exclude; n++) z*=z0;

    for(const_iterator ti=expansion.begin()+exclude;ti!=expansion.end();ti++) {
        B+=(*ti)*z;
        z*=z0;
    }

    return B0*B;
}

Vector3D MultipoleField::GetBFieldAt (const Point3D& x, double t) const
{
    Complex B = GetField2D(x.x,x.y);
    return Vector3D(B.imag(),B.real(),0);
}

Vector3D MultipoleField::GetEFieldAt (const Point3D& x, double t) const
{
    return Vector3D(0,0,0);
}

Vector3D MultipoleField::GetForceAt (const Point3D& x, const Vector3D& v, double q, double t) const
{
    Complex B = GetField2D(x.x,x.y);
    double qBx=q*B.imag();
    double qBy=q*B.real();
    return Vector3D(-qBy*v.z,qBx*v.z,qBy*v.x-qBx*v.y);
}

void MultipoleField::RotateY180 ()
{
    // To rotate the field 180 degrees, we perform the following
    //			bn		an
    // n-odd	-		+
    // n-even	+		-
    bool even=false;
    for(size_t n=0; n<expansion.size(); n++, even=!even) {
        Complex& z = expansion[n];
        if(even)
            z=Complex(z.real(),-z.imag());
        else
            z=Complex(-z.real(),z.imag());
    }
}

void MultipoleField::PrintField (std::ostream& os) const
{
    os<<"B0: "<<B0<<" Tesla\n";
    copy(expansion.begin(),expansion.end(),ostream_iterator<Complex>(os,"\n"));
    os<<endl;
}

void MultipoleField::SetComponent (int np, double bn, double an, double r0)
{
    assert(B0!=0);

    if(np+1>expansion.size())
        expansion.resize(np+1,Complex(0,0));

    expansion[np]=Complex(bn,an)/(B0*pow(r0,np));
}

Complex MultipoleField::GetComponent (int np, double r0) const
{
    return B0*GetCoefficient(np,r0);
}

Complex MultipoleField::GetCoefficient (int np, double r0) const
{
    if(np+1>expansion.size())
        const_cast<MultipoleField&>(*this).expansion.resize(np+1,Complex(0,0));

    return expansion[np]*pow(r0,np);
}

void MultipoleField::SetCoefficient (int np, const Complex& b, double r0)
{
    if(np+1>expansion.size())
        expansion.resize(np+1,Complex(0,0));
    expansion[np]=b/pow(r0,np);
}

