/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/11/23 10:16:11 $
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef BPM_h
#define BPM_h 1

#include "merlin_config.h"
#include <set>

// Monitor
#include "AcceleratorModel/StdComponent/Monitor.h"
// AMBufferManager
#include "AcceleratorModel/ActiveMonitors/AMBufferManager.h"
// Measurement
#include "NumericalUtils/Measurement.h"

class Bunch;
class BPM;
class ComponentTracker;

// Class representing a Beam Position Monitor BPM
// BPMs are active diagnostics which record the 
// following bunch properties:
//
//    x offset
//    y offset
//    total charge
//    ct
//

class BPM : public Monitor
{
public:

    // data structure for BPM data
    struct Data
    {
        double ct;
        Measurement x;
        Measurement y;
    };

    //	Buffer used by  BPM objects to record the results of a
    //	measurement. A single buffer can be associated with many
    //	BPM objects. When a BPM makes a measurement, it sends
    //	the results to its associated buffer using a BPM::Data
    //	struct via the Record() method.

    class Buffer
    {
    public:
        virtual ~Buffer ();
        virtual void Record (const BPM& aBPM, const Data& data) = 0;
    };

	// BPM::Response objects are used to modify the (x,y) measurement
	// obtained from the Bunch object. They are intended to mimic the
	// electrical (spacial) response of the monitor.

	class Response {
	public:
		virtual ~Response() {};
		virtual void Apply(Data*) =0;
	};

    typedef AMBufferManager<BPM,Buffer,Data> BufferManager;

    //	BPM constructing taking the id, the length and
    //	measurement point of the device.
    explicit BPM (const string& id, double len = 0, double mpos = 0);

    //	Set the resolution (rms noise levels) for the BPM.
    void SetResolution (double xr, double yr);
    double SetRes (double r) {
        SetResolution(r,r); return r;
    }

    //	Sets the scale factor for x and y planes (default scale
    //	=1).
    void SetScale (double xs, double ys);

	// Sets the response for the BPM. Returns the
	// original Response object (or NULL).
	Response* SetResponse(Response*);

    //	Measure the beam centroid of bunch.
    virtual void MakeMeasurement (const Bunch& aBunch);

    //	Returns the index for a BPM.
    virtual int GetIndex () const;

    //	Returns the type string "BPM".
    virtual const string& GetType () const;

    //	Primary tracking interface. Prepares the specified
    //	Tracker object for tracking this component.
    virtual void PrepareTracker (ComponentTracker& aTracker);

    //	Virtual constructor.
    virtual ModelElement* Copy () const;

    //	Sets the buffer to aBuffer. If a null pointer is passed,
    //	then no specific buffer is associated with this BPM, and
    //	the default buffer will be used (if not null).
    void AddBuffer (Buffer* abuffer);
    bool RemoveBuffer (Buffer* aBuffer);

    //	Removes all buffers associated with this BPM.
    void ClearAllBuffers ();

    //	Set the default buffer to be used by all BPMs, unless
    //	they are associated with a specific buffer.
    static void SetDefaultBuffer (Buffer* buffer);

    static const int ID;
    static bool generate_noise;

private:
    double res_x;
    double res_y;
    double scale_x;
    double scale_y;

	Response* itsResponse;

    BufferManager buffers;

    bool TakeData () const;
};

inline BPM::Buffer::~Buffer ()
{}

inline BPM::BPM (const string& id, double len, double mpos)
        : Monitor(id,len,mpos),res_x(0),res_y(0),
        scale_x(1),scale_y(1),itsResponse(0),buffers()
{}

inline void BPM::SetResolution (double xr, double yr)
{
    res_x=xr;
    res_y=yr;
}

inline void BPM::AddBuffer (BPM::Buffer* abuffer)
{
    buffers.AddBuffer(abuffer);
}

inline bool BPM::RemoveBuffer (BPM::Buffer* aBuffer)
{
    return buffers.RemoveBuffer(aBuffer);
}

inline void BPM::ClearAllBuffers ()
{
    buffers.ClearAllBuffers();
}

inline void BPM::SetDefaultBuffer (BPM::Buffer* buffer)
{
    BufferManager::SetDefaultBuffer(buffer);
}

inline bool BPM::TakeData () const
{
    return IsActive() && !buffers.empty();
}

inline BPM::Response* BPM::SetResponse(BPM::Response* r)
{
	Response* oldr = itsResponse;
	itsResponse = r;
	return oldr;
}


#endif
