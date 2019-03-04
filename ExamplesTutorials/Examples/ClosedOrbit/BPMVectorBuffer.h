/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "BPM.h"

class BPMVectorBuffer: public BPM::Buffer
{
public:
	BPMVectorBuffer()
	{
	}

	vector<BPM::Data> BPMReading;

	void Record(const BPM& aBPM, const BPM::Data& data)
	{
		BPMReading.push_back(data);
	}

	void Clear()
	{
		BPMReading.clear();
	}
};
