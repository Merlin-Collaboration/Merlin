/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (2000)
 * 
 * file Merlin\AcceleratorModel\StdComponent\SectorBend.cpp
 * last modified 04/04/01 15:24:08
 */

/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 *
 * Copyright (c) 2000 by The Merlin Collaboration.  
 * ALL RIGHTS RESERVED. 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */



// SectorBend
#include "AcceleratorModel/StdComponent/SectorBend.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"



#include "NumericalUtils/PhysicalConstants.h"
#include <utility>


// Class SectorBend::PoleFace




SectorBend::PoleFace::PoleFace (double angle, double f_int, double hg)
        : rot(angle),fint(f_int),hgap(hg)
{
}


// Class SectorBend::PoleFaceInfo





void SectorBend::PoleFaceInfo::Copy (const PoleFaceInfo& rhs)
{
    entrance = rhs.entrance ? new PoleFace(*(rhs.entrance)) : 0;
    if(rhs.entrance==rhs.exit)
        exit = entrance;
    else
        exit = rhs.exit ? new PoleFace(*(rhs.exit)) : 0;
}

// Class SectorBend

const int SectorBend::ID = UniqueIndex();


SectorBend::SectorBend (const string& id, double len, double h, double b0)
        : TAccCompGF<ArcGeometry,MultipoleField>(id,new ArcGeometry(len,h),new MultipoleField(0,b0,1.0,false))
{
}



int SectorBend::GetIndex () const
{
    return ID;
}

const string& SectorBend::GetType () const
{
    _TYPESTR(SectorBend);
}

void SectorBend::RotateY180 ()
{
    GetField().RotateY180();
    double h = GetGeometry().GetCurvature();
    GetGeometry().SetCurvature(-h);
    std::swap(pfInfo.entrance,pfInfo.exit);
}

void SectorBend::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

ModelElement* SectorBend::Copy () const
{
    return new SectorBend(*this);
}

void SectorBend::SetPoleFaceInfo (PoleFace* entr, PoleFace* exit)
{
    pfInfo.SetInfo(entr,exit);
}

void SectorBend::SetPoleFaceInfo (PoleFace* pf)
{
    pfInfo.SetInfo(pf);
}

double SectorBend::GetMatchedMomentum (double q) const
{
    using namespace PhysicalConstants;
    double h = GetGeometry().GetCurvature();
    double By = GetField().GetCoefficient(0).real();
    By*=GetField().GetFieldScale();
    return 1.0e-09*SpeedOfLight*By/h/q;
}

double SectorBend::GetB1 () const
{
    return GetField().GetComponent(1).real();
}

void SectorBend::SetB1 (double b1)
{
    GetField().SetComponent(1,b1);
}

