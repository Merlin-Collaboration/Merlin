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

// RandomNG
#include "Random/RandomNG.h"
// Bunch
#include "BeamModel/Bunch.h"
// RMSProfileMonitor
#include "AcceleratorModel/ActiveMonitors/RMSProfileMonitor.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

#define _ADDNOISE(var) if((var).error!=0) (var).value+=RandomNG::normal(0.0,(var).error*(var).error)

const int RMSProfileMonitor::ID = UniqueIndex();

int RMSProfileMonitor::GetIndex () const
{
    return ID;
}

const string& RMSProfileMonitor::GetType () const
{
    _TYPESTR(RMSProfileMonitor)
}

void RMSProfileMonitor::MakeMeasurement (const Bunch& aBunch)
{
    if(!buffers.empty() && IsActive()) {
        PSmoments2D profile;
        aBunch.GetProjectedMoments(ps_X,ps_Y,profile);

        Data mdat;

        mdat.x0.value = profile.mean(0);
        mdat.x0.error = res_x;
        mdat.y0.value = profile.mean(1);
        mdat.y0.error = res_y;
        mdat.xrms.value = profile.std(0);
        mdat.xrms.error = res_x;
        mdat.yrms.value = profile.std(1);
        mdat.yrms.error = res_y;

        // angled measurement
        if(uangle==0) {
            mdat.u0.value = mdat.x0.value;
            mdat.urms.value = mdat.xrms.value;
        }
        else if(uangle==pi/2) {
            mdat.u0.value = mdat.y0.value;
            mdat.urms.value = mdat.yrms.value;
        }
        else {
            double cosu = cos(uangle);
            double sinu = sin(uangle);
            mdat.u0.value = cosu*profile.mean(0)-sinu*profile.mean(1);
            mdat.urms.value = cosu*profile.std(0)-sinu*profile.std(1);
        }

        mdat.u0.error = res_u;
        mdat.urms.error = res_u;

        // Add noise

        _ADDNOISE(mdat.x0);
        _ADDNOISE(mdat.y0);
        _ADDNOISE(mdat.u0);
        _ADDNOISE(mdat.xrms);
        _ADDNOISE(mdat.yrms);
        _ADDNOISE(mdat.urms);

        buffers.SendToBuffers(*this,mdat);
    }
}

void RMSProfileMonitor::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Monitor)
}

ModelElement* RMSProfileMonitor::Copy () const
{
    return new RMSProfileMonitor(*this);
}

