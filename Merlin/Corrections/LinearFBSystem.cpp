/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
#include "TLAS/TLASimp.h"

// LinearFBSystem
#include "Corrections/LinearFBSystem.h"
using namespace TLAS;

namespace {

double ChannelRMS(const ROChannelArray& channels)
{
    RealVector v(channels.Size());
    channels.ReadAll(v);
    return sqrt(v*v/v.size());
}

};

LinearFBSystem::LinearFBSystem (std::vector<ROChannel*>& sigs, std::vector<RWChannel*>& acts, const RealMatrix& M)
        : gain(1.0),signals(sigs),actuators(acts),setpoints(0.0,sigs.size()),cached_actuators(0),Mi(0),
        actuatorQueue(0)
{
    SetResponseMatrix(M);
}

LinearFBSystem::LinearFBSystem (ROChannelArray& sigs, RWChannelArray& acts, const RealMatrix& M)
        : gain(1.0),signals(sigs),actuators(acts),setpoints(),cached_actuators(0),Mi(0),
        actuatorQueue(0)
{
    setpoints.redim(signals.Size());
    SetResponseMatrix(M);
}

LinearFBSystem::~LinearFBSystem ()
{
    if(Mi!=0) delete Mi;
    if(actuatorQueue!=0) delete actuatorQueue;
    if(cached_actuators!=0) delete cached_actuators;
}



void LinearFBSystem::SignalsToSetpoints ()
{
    signals.ReadAll(setpoints);
}

void LinearFBSystem::StoreActuators () const
{
    // Note that this stores the current actuator setting, regardless
    // of the state of actuatorQueue.
    if(cached_actuators==0)
        cached_actuators = new RealVector(GetNumActuators());
    actuators.ReadAll(*cached_actuators);
}

void LinearFBSystem::RestoreActuators ()
{
    assert(cached_actuators!=0);
    actuators.WriteAll(*cached_actuators);
}

void LinearFBSystem::SetResponseMatrix (const RealMatrix& M)
{
    assert(signals.Size()==M.nrows() && actuators.Size()==M.ncols());
    if(Mi!=0)
        delete Mi;
    Mi = new SVDMatrix<double>(M);
}

void LinearFBSystem::SetGain (double g)
{
    gain=g;
}

double LinearFBSystem::GetGain () const
{
    return gain;
}

void LinearFBSystem::Apply ()
{
    RealVector S(signals.Size());
    signals.ReadAll(S);
    S-=setpoints;

    RealVector newActs(actuators.Size());
    actuators.ReadAll(newActs);
    newActs += -gain*(*Mi)(S);

    if(actuatorQueue) {
        actuatorQueue->push(newActs);
        actuators.WriteAll(actuatorQueue->front());
        actuatorQueue->pop();
    }
    else
        actuators.WriteAll(newActs);
}

void LinearFBSystem::SetSetpoints (const RealVector& S0)
{
    assert(S0.size()==setpoints.size());
    setpoints=S0;
}

double LinearFBSystem::GetActuatorRMS () const
{
    return ChannelRMS(actuators);
}

double LinearFBSystem::GetSignalRMS () const
{
    RealVector S(signals.Size());
    signals.ReadAll(S);
    S-=setpoints;
    return sqrt(S*S/S.size());
}

int LinearFBSystem::GetNumSignals () const
{
    return signals.Size();
}

int LinearFBSystem::GetNumActuators () const
{
    return actuators.Size();
}

void LinearFBSystem::SetPulseDelay(int n)
{
    if(actuatorQueue!=0) {
        delete actuatorQueue;
        actuatorQueue = 0;
    }

    assert(n>0);

    if(n==1)
        return;

    actuatorQueue = new std::queue<RealVector>();

    RealVector zeros(0.0,actuators.Size());
    while((--n)>0)
        actuatorQueue->push(zeros);
}

