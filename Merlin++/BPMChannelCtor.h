/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BPMChannelCtor_h
#define BPMChannelCtor_h 1

#include "merlin_config.h"
#include "ChannelServer.h"
#include "BPMDataBufferServer.h"

class BPMChannel;

/**
 *	Specialised constructor for BPM channels. The type of
 *	channel constructed depends on the value of the
 *	inherited key string (='X' or 'Y').
 */

class BPMChannelCtor: public ChannelServer::ChannelCtor
{
public:

	/**
	 * Constructor taking either 'x' or 'y' to designate
	 * which plane of BPM channel to construct.
	 */
	BPMChannelCtor(char xy);

	/**
	 *	Constructs and returns a BPMChannel.
	 *	@return Constructed BPMChannel
	 */
	virtual ROChannel* ConstructRO(ModelElement* anElement);

	/**
	 *	Returns a nullptr.
	 *	@return nullptr
	 */
	virtual RWChannel* ConstructRW(ModelElement* anElement);

	static BPMDataBufferServer theServer;
};

#endif
