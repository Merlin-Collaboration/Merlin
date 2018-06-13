/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "BPMDataBufferServer.h"

BPMDataBuffer* BPMDataBufferServer::GetDataBuffer(BPM* bpm, bool create)
{
	map<BPM*, BPMDataBuffer>::iterator b = dataBuffers.find(bpm);

	if(b != dataBuffers.end())
	{
		return &(b->second);
	}

	if(!create)
	{
		return (BPMDataBuffer *) nullptr;
	}

	// make new buffer
	BPMDataBuffer* newbuf = &(dataBuffers[bpm]);
	newbuf->SetID(bpm->GetName());
	bpm->AddBuffer(newbuf);
	return newbuf;
}
