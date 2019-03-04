/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cassert>
#include "BPM.h"
#include "BPMChannelCtor.h"
#include "BPMChannel.h"

BPMDataBufferServer BPMChannelCtor::theServer;

BPMChannelCtor::BPMChannelCtor(char xy) :
	ChannelCtor("BPM", std::string(1, xy))
{
	assert(xy == 'X' || xy == 'Y');
}

ROChannel* BPMChannelCtor::ConstructRO(ModelElement* anElement)
{
#ifndef NDEBUG
	BPM* bpm = dynamic_cast<BPM*>(anElement);
	assert(bpm != 0);
#else
	BPM* bpm = static_cast<BPM*>(anElement);
#endif

	BPMDataBuffer* dataBuffer = theServer.GetDataBuffer(bpm);
	return new BPMChannel(key[0], dataBuffer);
}

RWChannel* BPMChannelCtor::ConstructRW(ModelElement* anElement)
{
	return (RWChannel *) nullptr;
}
