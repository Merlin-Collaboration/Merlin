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
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_RTMap
#define _h_RTMap 1

#include "BasicTransport/RMap.h"

// class RTMap
//
// A second-order TRANSPORT map for a PSvector
// class RTMap represents both the first-order R and second-order T TRANSPORT
// matrices. For efficiency, only non-zero terms are stored.

class RTMap : public RMap {
private:

class Tijk : public RMap::Rij {
    public:
        Tijk(int i1, int j1, int k1, double val =0)
                : Rij(i1,j1,val),k(k1) {}

        void Apply(const PSvector& orig,PSvector& res) const {
            res[i]+=val*orig[j]*orig[k];
        }

        int k;
    };

    typedef std::vector<Tijk> NonLinearTermList;
    typedef NonLinearTermList::iterator itor;
    typedef NonLinearTermList::const_iterator const_itor;

public:

    // Construction
    RTMap() : RMap(), tterms() {
        tterms.reserve(8);
    }

    // Construct linear map from matrix
    RTMap(const RealMatrix& R) : RMap(R), tterms() {
        tterms.reserve(8);
    }

    // Adding terms. Note that no effort
    // is  made to check that a term has
    // not been added more than once.
    double& operator()(int i, int j, int k) {
        tterms.push_back(Tijk(i-1,j-1,k-1));
        return tterms.back().val;
    }
    double& operator()(int i, int j) {
        return RMap::operator()(i,j);
    }

    // Operating on a PSvector
    PSvector& Apply(PSvector& p) const;

    // Output
    void Print(std::ostream&) const;

private:

    NonLinearTermList tterms;
};


#endif
