/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BPMDataBuffer_h
#define BPMDataBuffer_h 1

#include "merlin_config.h"
#include "BPM.h"

/**
 *	A BPM buffer which stores a single (x,y) pair from a
 *	BPM. The cached data is updated on each call to Record().
 *   A BPMbuffer can contain an additional offset.
 */

class BPMDataBuffer: public BPM::Buffer
{
public:

	/**
	 *	Constructor
	 */
	BPMDataBuffer();

	/**
	 * Record the data to the buffer.
	 */
	virtual void Record(const BPM& aBPM, const BPM::Data& data);

	/**
	 * Set the ID of the buffer.
	 */
	void SetID(const string& anID);

	/**
	 * Set buffer offsets (note: this is a temporary solution
	 * to setting arbitrary offsets to BPMs)
	 */
	void SetOffsets(double x, double y);

private:

	double x;
	double x_off;

	double y;
	double y_off;

	string id;

	friend class BPMChannel;
};

inline BPMDataBuffer::BPMDataBuffer() :
	x(0), x_off(0), y(0), y_off(0), id()
{
}

inline void BPMDataBuffer::SetID(const string& anID)
{
	id = anID;
}

inline void BPMDataBuffer::SetOffsets(double x, double y)
{
	x_off = x;
	y_off = y;
}

#endif
