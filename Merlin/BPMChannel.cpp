/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cassert>
#include "BPMDataBuffer.h"
#include "BPMChannel.h"

BPMChannel::BPMChannel(char XorY, BPMDataBuffer* dataBuff) :
	xy(XorY), itsData(dataBuff)
{
	assert((xy == 'X' || xy == 'Y') && (dataBuff != 0));
}

std::string BPMChannel::GetID() const
{
	return std::string("BPM." + itsData->id + '.' + xy);
}

double BPMChannel::Read() const
{
	if(xy == 'X')
	{
		return itsData->x + itsData->x_off;
	}
	else
	{
		return itsData->y + itsData->y_off;
	}
}

void BPMChannel::SetBPMOffset(double x)
{
	if(xy == 'X')
	{
		itsData->SetOffsets(x, 0);
	}
	else
	{
		itsData->SetOffsets(0, x);
	}
}
