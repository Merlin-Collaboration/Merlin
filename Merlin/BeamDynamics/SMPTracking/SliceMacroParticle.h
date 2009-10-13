/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_SliceMacroParticle
#define _h_SliceMacroParticle 1

#include "BeamModel/PSTypes.h"
#include "EuclideanGeometry/Space2D.h"
#include <iostream>
#include <cmath>

namespace SMPTracking {

typedef TPSMoments<2> PSmoments4D;

class SliceMacroParticle : public PSmoments4D {
public:

    explicit SliceMacroParticle(double q=0);
    SliceMacroParticle(const PSmoments& sigma, double ct, double dp, double q);

    // Macroparticle charge
    double Q() const { return q; }


    // weighted centroid values
    double GetChargeWeightedCentroid(PScoord i) const {
        return q*mean(i);
    }
    Point2D GetChargeWeightedCentroid(PScoord i, PScoord j) const {
        return Point2D(q*mean(i),q*mean(j));
    }
    PSvector GetChargeWeightedCentroid() const {
        PSvector x(*this);
        x*=q;
        return x;
    }

    PSvector& GetCentroid() { return *this; }
    const PSvector& GetCentroid() const { return *this; }

    // IO methods
    void Read(std::istream&);
    void Write(std::ostream&) const;

    friend bool operator<(const SliceMacroParticle& p1, const SliceMacroParticle& p2) {
        return p1.ct()<p2.ct();
    }

private:
    double q;
};

inline std::ostream& operator<<(std::ostream& os, const SliceMacroParticle& p)
{
    p.Write(os); return os;
}

inline std::istream& operator>>(std::istream& is, SliceMacroParticle& p)
{
    p.Read(is); return is;
}

}; // end namespace SMPTracking

#endif
