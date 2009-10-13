/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/Miscellaneous/CorrectorWinding.h"

CorrectorWinding::CorrectorWinding (RectMultipole& aMagnet)
        : ModelElement(aMagnet.GetName()),magnet(&aMagnet)
{}

void CorrectorWinding::SetBx (double value)
{
    if((*magnet).GetField().GetFieldScale()==0)
        return;

    double l=magnet->GetLength();
    (*magnet).GetField().SetComponent(0,0,value/l);
}

void CorrectorWinding::SetBy (double value)
{
    if((*magnet).GetField().GetFieldScale()==0)
        return;

    double l=magnet->GetLength();
    (*magnet).GetField().SetComponent(0,value/l,0);
}

double CorrectorWinding::GetBx () const
{
    return (*magnet).GetField().GetComponent(0).imag()*magnet->GetLength();
}

double CorrectorWinding::GetBy () const
{
    return (*magnet).GetField().GetComponent(0).real()*magnet->GetLength();
}

const string& CorrectorWinding::GetType () const
{
    _TYPESTR(CorrectorWinding)
}

ModelElement* CorrectorWinding::Copy () const
{
    return new CorrectorWinding(*this);
}

