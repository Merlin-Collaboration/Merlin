/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BPMDataBufferServer_h
#define BPMDataBufferServer_h 1

#include "merlin_config.h"
#include <map>
#include "BPMDataBuffer.h"

/**
 *	Singleton instance of a storage class for BPMChannel
 *	Data objects. BPMChannelData are indexed by their
 *	associated BPM pointer.
 */

class BPMDataBufferServer
{
public:

	BPMDataBuffer* GetDataBuffer(BPM* bpm, bool create = true);
	void Dump(ostream& os);

private:

	std::map<BPM*, BPMDataBuffer> dataBuffers;
};

#endif
