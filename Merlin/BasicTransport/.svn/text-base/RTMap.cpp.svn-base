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

#include "BasicTransport/RTMap.h"
#include "NumericalUtils/MatrixPrinter.h"
#include "TLAS/LinearAlgebra.h"

void RTMap::Print(std::ostream& os) const
{
    // Print linear matrix first
    os<<"R:"<<endl;
    RMap::MatrixForm(os);

    os<<"\nT:\n";
    for(int n=0;n<6;n++) {
        RealMatrix M(6,6,0.0);
        for(const_itor t = tterms.begin(); t!=tterms.end(); t++) {
            if((t->i)==n) {
                if(t->j != t->k)
                    M(t->k,t->j) = M(t->j,t->k) = t->val/2;
                else
                    M(t->j,t->k) = t->val;
            }
        }
        ::MatrixForm(M,os);
        os<<endl;
    }
}

PSvector& RTMap::Apply(PSvector& X) const
{
    PSvector Y(0);

    // linear map
    RMap::Apply(X,Y);

    // non-linear map
    for(const_itor t = tterms.begin(); t!=tterms.end(); t++)
        t->Apply(X,Y);

    return X=Y;
};

