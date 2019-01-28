/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "merlin_config.h"
#include "Bunch.h"
#include "BPM.h"
#include "ComponentTracker.h"
#include "RandomNG.h"

#include "utils.h"

using namespace std;

bool BPM::generate_noise = true;
const int BPM::ID = UniqueIndex();

void BPM::SetScale(double xs, double ys)
{
	scale_x = xs;
	scale_y = ys;
}

void BPM::MakeMeasurement(const Bunch& aBunch)
{
	if(TakeData())
	{
		Point2D x0 = aBunch.GetProjectedCentroid(ps_X, ps_Y);

		if(!fequal(generate_noise, 0.0) && !fequal(res_x, 0.0))
		{
			x0.x += RandomNG::normal(0.0, res_x * res_x);
		}
		if(!fequal(generate_noise, 0.0) && !fequal(res_y, 0.0))
		{
			x0.y += RandomNG::normal(0.0, res_y * res_y);
		}

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
		{
			itsResponse->Apply(&mdat);
		}

		buffers.SendToBuffers(*this, mdat);
	}
}

int BPM::GetIndex() const
{
	return ID;
}

const string& BPM::GetType() const
{
	_TYPESTR(BPM)
}

void BPM::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, Monitor)
}

ModelElement* BPM::Copy() const
{
	return new BPM(*this);
}
