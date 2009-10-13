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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef PSmoments_h
#define PSmoments_h 1

#include "merlin_config.h"
#include "BeamModel/PSvector.h"
#include "NumericalUtils/TCovMtrx.h"
#include <iostream>
#include <iomanip>

typedef TCovMtrx<double,6> SigmaMatrix;
typedef TCovMtrx<double,4> SigmaMatrix4D;
typedef TCovMtrx<double,2> SigmaMatrix2D;

// template class that stores the first- and second-order
// phase space moments. The first-order moments are always
// stored as a six-vector (PSvector), while the template
// paramter defines the number of degrees of freedom for
// the second-order moments.

template<int N>
class TPSMoments : public PSvector, public TCovMtrx<double,2*N> {
public:

    typedef TCovMtrx<double,2*N> SigMtrx;

    TPSMoments()  :PSvector(0),SigMtrx() {}

    // mean added for backwards compatability
    // functions delegated to PSvector and SigMtrx
    double mean(int i) const {return operator[](i);}
    double& mean(int i) {return operator[](i);}

    // print the moments as a table (same as TRANSPORT output)
    void printFormatted(std::ostream& os, bool norm =true) const;
};

// template specialisation for 1 degree of freedom.
template<>
class TPSMoments<1> : public SigmaMatrix2D {
public:

    TPSMoments<1>() : SigmaMatrix2D() {
        data[0]=data[1]=0;
    }

double operator[](int i) const { return i<2? data[i] : 0;}
    double& operator[](int i) {return data[i];}

    // mean added for backwards compatability
    double mean(int i) const {
        return operator[](i);
    }

private:
    double data[2];
};

typedef TPSMoments<1> PSmoments2D;
typedef TPSMoments<3> PSmoments;
typedef std::vector<PSmoments> PSmomentsArray;

template <int N>
void TPSMoments<N>::printFormatted (std::ostream& os, bool normalised) const
{
    using std::setw;

    int oldp=os.precision(3);
    ios_base::fmtflags oflg=os.setf(ios::scientific,ios::floatfield);

    for(int i=0; i<2*N; i++) {
        os<<setw(12)<<right<<scientific<<this->mean(i);
        os<<setw(12)<<right<<scientific<<this->std(i);
        for(int j=1; j<2*N; j++) {
            if(j<=i) {
                if(normalised)
                    os<<setw(8)<<" ";
                else
                    os<<setw(12)<<" ";
            }
            else {
                os<<right;
                if(normalised)
                    os<<setw(8)<<fixed<<this->r_ij(i,j);
                else
                    os<<setw(12)<<scientific<<this->sig(i,j);
            }
        }
        os<<endl;
    }
    os.precision(oldp);
    os.flags(oflg);
}

#endif
