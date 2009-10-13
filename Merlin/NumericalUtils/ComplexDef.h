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

#ifndef _h_ComplexDef
#define _h_ComplexDef 1

#include <cmath>
#include <iostream>

class Complex {
    double _re,_im;
public:
    typedef double value_type;
    double real() const {return _re;}
    double imag() const {return _im;}
    void real(double re) { _re=re;}
    void imag(double im) { _im=im;}
    Complex(double re = 0, double im = 0):_re(re),_im(im){}
    Complex& operator+=(const Complex& rhs) {
        _re+=rhs._re;_im+=rhs._im;
        return *this;
    }
    Complex& operator-=(const Complex& rhs) {
        _re-=rhs._re;_im-=rhs._im;
        return *this;
    }
    Complex& operator*=(const Complex& rhs) {
        double x = _re*rhs._re-(_im)*rhs._im;
        double y = _re*rhs._im+(_im)*rhs._re;
        _re=x;_im=y;
        return *this;
    }
    Complex& operator/=(const Complex& rhs) {
        double d=rhs._re*rhs._re+rhs._im*rhs._im;
        double x = (_re*rhs._re+(_im)*rhs._im)/d;
        double y = (_re*rhs._im-(_im)*rhs._re)/d;
        _re=x;_im=y;
        return *this;
    }
    Complex& operator=(double rhs) {
        _re =rhs; _im=double(0); return *this;
    }
    Complex& operator+=(double rhs) {
        _re+=rhs; return *this;
    }
    Complex& operator-=(double rhs) {
        _re-=rhs; return *this;
    }
    Complex& operator*=(double rhs) {
        _re*=rhs;_im*=rhs; return *this;
    }
    Complex& operator/=(double rhs) {
        _re/=rhs;_im/=rhs; return *this;
    }
    friend Complex operator+(const Complex& lhs, const Complex& rhs) {
        return Complex(lhs._re+rhs._re,lhs._im+rhs._im);
    }
    friend Complex operator+(const Complex& lhs, double rhs) {
        return Complex(lhs._re+rhs,lhs._im);
    }
    friend Complex operator+(double lhs, const Complex& rhs) {
        return Complex(lhs+rhs._re,rhs._im);
    }
    friend Complex operator-(const Complex& lhs, const Complex& rhs) {
        return Complex(lhs._re-rhs._re,lhs._im-rhs._im);
    }
    friend Complex operator-(const Complex& lhs, double rhs) {
        return Complex(lhs._re-rhs,lhs._im);
    }
    friend Complex operator-(double lhs, const Complex& rhs) {
        return Complex(lhs-rhs._re,-rhs._im);
    }
    friend Complex operator*(const Complex& lhs, double rhs) {
        return Complex(rhs*lhs._re,rhs*lhs._im);
    }
    friend Complex operator*(double lhs, const Complex& rhs) {
        return operator*(rhs,lhs);
    }
    friend Complex operator*(const Complex& rhs, const Complex& lhs) {
        Complex t(rhs); return t*=lhs;
    }
    friend Complex operator/(const Complex& lhs, double rhs) {
        return Complex(lhs._re/rhs,lhs._im/rhs);
    }
    friend Complex operator/(double lhs, const Complex& rhs) {
        Complex tmp(lhs); return tmp/=rhs;
    }
    friend bool operator==(const Complex& lhs, const Complex& rhs) {
        return lhs._re==rhs._re && lhs._im==rhs._im;
    }
    friend bool operator!=(const Complex& lhs, const Complex& rhs) {
        return lhs._re!=rhs._re || lhs._im!=rhs._im;
    }
    friend bool operator==(const Complex& lhs, double rhs) {
        return lhs._re==rhs && lhs._im==double(0);
    }
    friend bool operator==(double lhs, const Complex& rhs) {
        return lhs==rhs._re && rhs._im==double(0);
    }
    friend bool operator!=(const Complex& lhs, double rhs) {
        return lhs._re!=rhs || lhs._im!=double(0);
    }

    friend bool operator!=(double lhs, const Complex& rhs) {
        return lhs!=rhs._re || rhs._im==double(0);
    }

    friend double abs(const Complex& rhs) {
        return ::sqrt((rhs._re)*(rhs._re)+(rhs._im)*(rhs._im));
    }
    friend double arg(const Complex& rhs) {
        return atan2(rhs._im,rhs._re);
    }
    friend bool operator<(const Complex&,const Complex&) {
        return false;
    }
    friend std::ostream& operator<<(std::ostream& os, const Complex& z) {
        return os<<z.real()<<"+"<<z.imag()<<"i";
    }
    friend Complex log(const Complex& z) {
        return Complex(::log(abs(z)),arg(z));
    }
};

#endif

