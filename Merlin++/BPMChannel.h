/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BPMChannel_h
#define BPMChannel_h 1

#include "merlin_config.h"
#include "Channels.h"

class BPMDataBuffer;

class BPMChannel: public ROChannel
{
public:

	/**
	 *	Constructor taking the BPMDataBuffer, and the type of
	 *	channel ('X' or 'Y').
	 */
	BPMChannel(char XorY, BPMDataBuffer* dataBuff);

	/**
	 *	Returns the ID of the channel (parameter).
	 *	@return Channel ID
	 */
	virtual std::string GetID() const;

	/**
	 *	Returns the current value of the parameter/attribute
	 *	associated with the channel.
	 *	@return Current value of the associated parameter/attribute
	 */
	virtual double Read() const;

	/**
	 * Special function for setting arbitrary offset
	 */
	void SetBPMOffset(double v);

private:

	char xy;
	BPMDataBuffer* itsData;
};

#endif
