/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cassert>
#include "TLASimp.h"

#include "LinearFBSystem.h"
using namespace TLAS;

namespace
{

double ChannelRMS(const ROChannelArray& channels)
{
	RealVector v(channels.Size());
	channels.ReadAll(v);
	return sqrt(v * v / v.size());
}

}

LinearFBSystem::LinearFBSystem(std::vector<ROChannel*>& sigs, std::vector<RWChannel*>& acts, const RealMatrix& M) :
	gain(1.0), signals(sigs), actuators(acts), setpoints(0.0, sigs.size()), cached_actuators(nullptr), Mi(nullptr),
	actuatorQueue(nullptr)
{
	SetResponseMatrix(M);
}

LinearFBSystem::LinearFBSystem(ROChannelArray& sigs, RWChannelArray& acts, const RealMatrix& M) :
	gain(1.0), signals(sigs), actuators(acts), setpoints(), cached_actuators(nullptr), Mi(nullptr), actuatorQueue(
		nullptr)
{
	setpoints.redim(signals.Size());
	SetResponseMatrix(M);
}

LinearFBSystem::~LinearFBSystem()
{
	if(Mi != nullptr)
	{
		delete Mi;
	}
	if(actuatorQueue != nullptr)
	{
		delete actuatorQueue;
	}
	if(cached_actuators != nullptr)
	{
		delete cached_actuators;
	}
}

void LinearFBSystem::SignalsToSetpoints()
{
	signals.ReadAll(setpoints);
}

void LinearFBSystem::StoreActuators() const
{
	// Note that this stores the current actuator setting, regardless
	// of the state of actuatorQueue.
	if(cached_actuators == nullptr)
	{
		cached_actuators = new RealVector(GetNumActuators());
	}
	actuators.ReadAll(*cached_actuators);
}

void LinearFBSystem::RestoreActuators()
{
	assert(cached_actuators != nullptr);
	actuators.WriteAll(*cached_actuators);
}

void LinearFBSystem::SetResponseMatrix(const RealMatrix& M)
{
	assert(signals.Size() == M.nrows() && actuators.Size() == M.ncols());
	if(Mi != nullptr)
	{
		delete Mi;
	}
	Mi = new SVDMatrix<double>(M);
}

void LinearFBSystem::SetGain(double g)
{
	gain = g;
}

double LinearFBSystem::GetGain() const
{
	return gain;
}

void LinearFBSystem::Apply()
{
	RealVector S(signals.Size());
	signals.ReadAll(S);
	S -= setpoints;

	RealVector newActs(actuators.Size());
	actuators.ReadAll(newActs);
	newActs += -gain * (*Mi)(S);

	if(actuatorQueue)
	{
		actuatorQueue->push(newActs);
		actuators.WriteAll(actuatorQueue->front());
		actuatorQueue->pop();
	}
	else
	{
		actuators.WriteAll(newActs);
	}
}

void LinearFBSystem::SetSetpoints(const RealVector& S0)
{
	assert(S0.size() == setpoints.size());
	setpoints = S0;
}

double LinearFBSystem::GetActuatorRMS() const
{
	return ChannelRMS(actuators);
}

double LinearFBSystem::GetSignalRMS() const
{
	RealVector S(signals.Size());
	signals.ReadAll(S);
	S -= setpoints;
	return sqrt(S * S / S.size());
}

int LinearFBSystem::GetNumSignals() const
{
	return signals.Size();
}

int LinearFBSystem::GetNumActuators() const
{
	return actuators.Size();
}

void LinearFBSystem::SetPulseDelay(int n)
{
	if(actuatorQueue != nullptr)
	{
		delete actuatorQueue;
		actuatorQueue = nullptr;
	}

	assert(n > 0);

	if(n == 1)
	{
		return;
	}

	actuatorQueue = new std::queue<RealVector>();

	RealVector zeros(0.0, actuators.Size());
	while((--n) > 0)
	{
		actuatorQueue->push(zeros);
	}
}
