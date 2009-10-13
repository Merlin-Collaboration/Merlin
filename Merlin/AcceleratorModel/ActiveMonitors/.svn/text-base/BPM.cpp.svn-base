/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/11/14 09:57:11 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#include "merlin_config.h"
// Bunch
#include "BeamModel/Bunch.h"
// BPM
#include "AcceleratorModel/ActiveMonitors/BPM.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// RandomNG
#include "Random/RandomNG.h"

using namespace std;

bool BPM::generate_noise = true;
const int BPM::ID = UniqueIndex();

void BPM::SetScale (double xs, double ys)
{
    scale_x = xs;
    scale_y = ys;
}

void BPM::MakeMeasurement (const Bunch& aBunch)
{
    if(TakeData()) {
        Point2D x0 = aBunch.GetProjectedCentroid(ps_X,ps_Y);

        if(generate_noise && res_x!=0)
            x0.x += RandomNG::normal(0.0,res_x*res_x);
        if(generate_noise && res_y!=0)
            x0.y += RandomNG::normal(0.0,res_y*res_y);

        x0.x *= scale_x;
        x0.y *= scale_y;

        Data mdat;
        mdat.x.value = x0.x;
        mdat.x.error = res_x;
        mdat.y.value = x0.y;
        mdat.y.error = res_y;
        //		mdat.q.value = aBunch.GetTotalCharge();
        //		mdat.q.error = res_q;
        mdat.ct = aBunch.GetReferenceTime();

		if(itsResponse)
			itsResponse->Apply(&mdat);

        buffers.SendToBuffers(*this,mdat);
    }
}

int BPM::GetIndex () const
{
    return ID;
}

const string& BPM::GetType () const
{
    _TYPESTR(BPM)
}

void BPM::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,Monitor)
}

ModelElement* BPM::Copy () const
{
    return new BPM(*this);
}

