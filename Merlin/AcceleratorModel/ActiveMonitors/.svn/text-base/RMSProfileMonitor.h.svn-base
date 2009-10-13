/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef RMSProfileMonitor_h
#define RMSProfileMonitor_h 1

#include "merlin_config.h"

#include "NumericalUtils/NumericalConstants.h"
#include <set>
// AMBufferManager
#include "AcceleratorModel/ActiveMonitors/AMBufferManager.h"
// Monitor
#include "AcceleratorModel/StdComponent/Monitor.h"
// Measurement
#include "NumericalUtils/Measurement.h"

class Bunch;
class ComponentTracker;

//	An RMS profile monitor. A profile monitor which mimicks
//	the action of a "wire scanner". A profile monitor can
//	measure the rms beam projection onto three planes:
//	horizontal (x), vertical (y) and a third so-called
//	u-wire, which is at a specifid angle to the x plane.

class RMSProfileMonitor : public Monitor
{
public:

    // Data structure for monitor data
    struct Data
    {
        Measurement x0;
        Measurement y0;
        Measurement u0;
        Measurement xrms;
        Measurement yrms;
        Measurement urms;
    };

    class Buffer
    {
    public:
        virtual ~Buffer ()
        {}
        virtual void Record (const RMSProfileMonitor& mon, const Data& dat) = 0;
    };

    typedef AMBufferManager< RMSProfileMonitor,Buffer,Data > BufferManager;


    RMSProfileMonitor (const string& id, double uphi = pi/4, double len = 0, double mpt = 0)
            : Monitor(id,len),res_x(0),res_y(0),res_u(0),uangle(uphi)
    {}

    void SetResolution (double rx, double ry, double ru)
    {
        res_x = rx;
        res_y = ry;
        res_u = ru;
    }

    void AddBuffer (Buffer* buffer)
    {
        buffers.AddBuffer(buffer);
    }

    bool RemoveBuffer (Buffer* aBuffer)
    {
        return buffers.RemoveBuffer(aBuffer);
    }

    void RemoveAllBuffers ()
    {
        buffers.ClearAllBuffers();
    }

    static void SetDefaultBuffer (Buffer* buffer)
    {
        BufferManager::SetDefaultBuffer(buffer);
    }

    //	Returns the unique index for this class of accelerator
    //	components.
    virtual int GetIndex () const;

    //	Returns the type string for this component.
    virtual const string& GetType () const;

    //	Pure virtual function. Makes a measurement on the
    //	supplied Beam object. Concrete diagnostics must supply
    //	this function.
    virtual void MakeMeasurement (const Bunch& aBunch);

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    static const int ID;

private:

    double res_x;
    double res_y;
    double res_u;
    double uangle;

    BufferManager buffers;
};

#endif
